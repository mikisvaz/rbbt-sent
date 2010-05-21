require 'sent'
require 'rbbt/sources/pubmed'
require 'rbbt/util/misc'
require 'rbbt/util/open'
require 'rbbt/bow/bow'
require 'rbbt/bow/dictionary'

require 'soap/wsdlDriver'
require 'stemmer'
require 'open4'
require 'progress-monitor'
require 'yaml'
require 'base64'

# Produce Stem lists, used by the Web Service
class String

  alias old_stem stem

  def self.reset_stem_list
    @@stem_dictionary = Hash.new
  end 

  # Extends the stem functionality so that is generates a dictionary of
  # stems. For each stem a list of words that reduce to it.
  def stem
    res = old_stem
    @@stem_dictionary[res] ||= Hash.new
    @@stem_dictionary[res][self] ||= 1
    res
  end

  # Returns the dictionary of recorded stems.
  def self.stem_list(dictionary = nil)
    stem_list = Hash.new
    @@stem_dictionary.each{|k,l|
      next if dictionary && !dictionary.include?(k)
      stem_list[k] = l.keys
    }
    stem_list
  end

  reset_stem_list
end

def try3times
  try = 0
  begin
    yield
  rescue
    try += 1
    retry if try <= 3
  end
end


module Sent

  class NoGenesError < StandardError; end
  class ProcessAbortedError < StandardError; end
  class WSError < StandardError; end

  def self.run_R(command)
    pid, stdin, stdout, stderr = Open4::popen4 "R --vanilla --slave"
    stdin.write "source('#{File.join(Sent.rdir,'matrix.R')}');\n"
    stdin.write "#{ command };\n"   
    stdin.close

    Process.wait pid

    raise Sent::ProcessAbortedError, "Error in R process: #{stdout.read + stderr.read}" if $?.exitstatus != 0
    result = stdout.read + stderr.read
    stdout.close
    stderr.close

    puts result if result != ""
    result
  end

  def self.mentions(org)
    ner   = Organism.ner(org, :rner)
    norm  = Organism.norm(org)
    pmids = Organism.literature(org)

    chunks = pmids.chunk(100)

    pmids = {}
    Progress.monitor("Finding gene-article associations in text", :step => 1000)
    chunks.each{|chunk|
      try3times do
        PubMed.get_article(chunk).each{|pmid, article|
          text = article.text

          mentions = ner.extract(text)

          Progress.monitor("Resolving mentions", :step => 1000)
          codes = mentions.collect{|mention| 
            matches = norm.match(mention)
            norm.select(matches,mention,text)
          }.flatten.uniq.sort

          codes.each{|code|
            pmids[code] ||= []
            pmids[code] << pmid
          }
        }
      end
    }
    pmids
  end

  def self.dictionary(associations, options = {})
    dict_options = {:low => 0.001, :hi => 0.65, :limit => 3000}.merge options
    dict = Dictionary::TF_IDF.new

    String.reset_stem_list

    associations.each do |gene, pmids|
      try3times do
        text = PubMed.get_article(pmids).collect{|pmid, article| article.text }.join("\n")
        dict.add(BagOfWords.count(text.bigrams))
      end
    end

    term_weigths = dict.weights(dict_options)
    terms        = term_weigths.keys.sort

    stems        = String.stem_list(terms.collect{|p| p.split(/ /) }.flatten.uniq)

    [term_weigths, stems]
  end
 

  def self.metadoc(associations, terms)
    "%s\t%s\n" % ["Entity", terms.sort * "\t"] +
    associations.collect do |gene, pmids|
      try3times do
        text = PubMed.get_article(pmids).collect{|pmid, article| article.text }.join("\n")
        "%s\t%s" % [gene, BagOfWords.features(text, terms.sort).join("\t")]
      end
    end * "\n"
  end

  def self.matrix(metadoc_file, genes)

    matrix = ""
    File.open(metadoc_file) do |f|
      matrix += f.gets
      f.read.each_line do |line|
        gene = line.match(/^(.*?)\t/)[1]
        matrix += line if genes.include? gene
      end
    end

    matrix
  end

  @@bionmf_wsdl = "http://bionmf.dacya.ucm.es/WebService/BioNMFWS.wsdl"   
  def self.NMF(matrix, out, k, executions = 10)
    driver = SOAP::WSDLDriverFactory.new( @@bionmf_wsdl).create_rpc_driver

    # Upload matrix
    nmf_matrix = driver.upload_matrix(
       matrix,   # matrix
       false,    # binary 
       true,     # column labels 
       true,     # row labels
       true,     # transpose
       "No",     # positive
       "No",     # normalization 
       "matrix") # Suggested name

    # Send several executions in parallel
    while !driver.done(nmf_matrix)
      sleep(5)
    end

    if driver.error(nmf_matrix) 
      error = driver.messages(nmf_matrix).join("\n")
      raise "Error pre-processing matrix!"  + driver.messages(nmf_matrix).join("\n")
    end

    threads = []
    error = nil
    executions.times{|i|
      threads << Thread.new(i){ |num|
        times = 3
        begin
          
          job_id = driver.standardNMF(
            nmf_matrix, # Matrix job
            "Standard", # Algorithm
            k,          # Factor Start
            k,          # Factor End
            1,          # Runs
            2000,       # Iterations
            40,         # Stop criteria
            0,          # Not used (nsnmf smoothness)
            false,      # extra info
            '')         # Suggested name

          while !driver.done(job_id)
            sleep(5)
          end

          if driver.error(job_id) 
            error = driver.messages(job_id).join("\n")
            raise "Error in NMF"  + driver.messages(job_id).join("\n")
          end

          results =  driver.results(job_id)
          
          File.open(out + ".matrix_w.#{num}",'w') do |f| 
            f.write Base64.decode64 driver.result(results[0]) #.sub(/\t(.*)\t$/,'\1')
          end

          File.open(out + ".matrix_h.#{num}",'w') do |f|
            f.write Base64.decode64 driver.result(results[1]) #.sub(/\t(.*)\t$/,'\1')
          end

          driver.clean(job_id)
        rescue Sent::ProcessAbortedError
          puts "Process aborted for #{ num }"
          driver.abort(job_id)
        rescue Timeout::Error
          if times > 0
            times -= 1
            sleep 2
            retry
          else
            raise Sent::ProcessAbortedError, "NMF Execution #{ num } timed out"
          end
        rescue Exception
          puts $!.message
          if times > 0
            times -= 1
            puts "Retrying thread #{ num }"
            retry
          else

            puts "NMF Execution #{ num } Produced Exception"
            puts $!.class
            puts $!.message
            puts $!.backtrace
            raise Sent::ProcessAbortedError, "NMF Execution #{ num } Produced Exception"
          end
        ensure
          Thread.exit
        end
      }
      sleep 1

    }

    # Allow threads to be aborted
    aborted = false
    old_int = Signal.trap("INT") do
      STDERR.puts "Killing threads"
      threads.each{|t| t.raise Sent::ProcessAbortedError, "Process Aborted"}
      aborted = true
    end

    threads.each { |aThread|  aThread.join }

    Signal.trap("INT", old_int)
    driver.clean(nmf_matrix)

    if aborted
      raise Sent::ProcessAbortedError, "Process Aborted"
    end

    if error
      raise Exception, "Error in NMF:\n" + error
    end
    
    run_R("SENT.join.results('#{ out }')")

    FileUtils.rm Dir.glob(out + '.matrix_*.*')
  end

  def self.CCC(matrix, kstart, kend)
    raise "Error in range: #{ kstart } to #{ kend }" if kstart >= kend

    driver = SOAP::WSDLDriverFactory.new( @@bionmf_wsdl).create_rpc_driver

    # Prepare matrix for processing
    nmf_matrix = driver.upload_matrix(File.open(matrix).read)
    driver.preprocess(nmf_matrix,1,"No","No", true, true)    

    job_id = driver.sample_classification(nmf_matrix,kstart.to_i,kend.to_i,10)

    aborted = false
    old_int = Signal.trap("INT") do
      puts "Aborting bestK process"
      driver.abort(job_id)
      aborted = true
    end

    while (status = driver.status(job_id)) == 0
      sleep(5)
    end

    driver.clean_matrix(nmf_matrix)
    Signal.trap("INT", old_int)

    if aborted
      raise Sent::ProcessAbortedError, "Process Aborted"
    end

    if status == -1 
      raise Sent::WSError, "Error processing matrix:\n" + driver.info(job_id)
    end

    results = driver.results(job_id)
    text = driver.get_result(results[0])
    text.split(/\n/s).last.split(/\t/)
  end

  def self.analyze(prefix, output, clusters = nil, num_words = 15)

    FileUtils.rm Dir.glob(output + '*.words') + Dir.glob(output + '*.genes')
    run_R("SENT.analyze('#{ prefix }', '#{ output }', #{clusters || 'NULL'}, #{num_words})")
    words = Dir.glob(output + '*.words').sort.collect{|f| Open.read(f).split(/\n/)}
    genes = Dir.glob(output + '*.genes').sort.collect{|f| Open.read(f).split(/\n/)}

    groups = []
    words.zip(genes).each{|p|
      groups << {:words => p[0], :genes => p[1]}
    }

    groups
  end

  def self.literature_index(pmids, outfile)

    index = Ferret::Index::Index.new(:path => outfile)

    index.field_infos.add_field(:title,  :index => :yes, :boost => 0.67)
    index.field_infos.add_field(:abstract, :index => :yes, :boost => 0.33)

    Progress.monitor("Building index for #{pmids.length} articles")
    pmids.each{|pmid|
      begin
        article = PubMed.get_article(pmid)
        abstract = article.abstract
        title    = article.title

        abstract_content = BagOfWords.terms(abstract).collect{|w,n| (1..n).collect{ w }}.flatten.join(" ")
        title_content    = BagOfWords.terms(title).collect{|w,n| (1..n).collect{ w }}.flatten.join(" ")

        index << {:id => pmid, :abstract => abstract_content, :name => title_content}
      rescue Exception
        puts $!.backtrace
        puts $!.message
      end

    }
    index.close
  end
  
  def self.search_index(words, index)
    index = Ferret::Index::Index.new(:path => index)

    ranks = []
    index.search_each("#{ words.collect{|w| w.stem}.join(" ") }", :limit => 8000) do |id,score|
      next unless score > 0.0001
      ranks << [index[id][:id],score]
    end

    ranks
  end
end
