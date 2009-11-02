$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require 'yaml'
require 'fileutils'
require 'stemmer'
require 'ferret'

module Sent

  class NoConfig < Exception; end

  @@rootdir = File.dirname(File.dirname(__FILE__))
  @@datadir = @@workdir = @@tmpdir = nil

  def self.load_config
    if File.exist?File.join(@@rootdir, 'sent.config')
      config = YAML.load_file(File.join(@@rotdir, 'sent.config'))
      if config.is_a? Hash
        @@datadir  = config['datadir'] if config['datadir'] 
        @@workdir = config['workdir'] if config['workdir']
        @@tmpdir   = config['tmpdir'] if config['tmpdir']
      end
    end



    if File.exist?(File.join(ENV['HOME'], '.sent'))
      config = YAML.load_file(File.join(ENV['HOME'], '.sent') )
      if config.is_a? Hash
        @@datadir  = config['datadir'] if config['datadir'] 
        @@workdir = config['workdir'] if config['workdir']
        @@tmpdir   = config['tmpdir'] if config['tmpdir']
      end
    end

    if @@datadir.nil?  || @@workdir.nil? || @@tmpdir.nil?
      raise Sent::NoConfig, "sent not configured. Edit #{File.join(@@rootdir, 'sent.config')} or $HOME/.sent"
    end


    FileUtils.mkdir_p @@datadir  unless File.exist? @@datadir
    FileUtils.mkdir_p @@workdir unless File.exist? @@workdir
    FileUtils.mkdir_p @@tmpdir   unless File.exist? @@tmpdir

  end

  def self.datadir
    @@datadir
  end 
  def self.workdir
    @@workdir
  end 
  def self.tmpdir
    @@tmpdir
  end 

  def self.rootdir
    @@rootdir
  end


  def self.rdir
    File.join(@@rootdir, 'R')
  end

  self.load_config
end


require 'rbbt/sources/pubmed'
require 'rbbt/util/misc'
require 'rbbt/util/open'

require 'soap/wsdlDriver'



require 'stemmer'
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
    @@stem_dictionary[res][self] = 1
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


module Sent

  class NoGenes < Exception; end
  class ProcessAborted < Exception; end
  class WSError < Exception; end

  require 'open4'
  def self.run_R(command)
    pid, stdin, stdout, stderr = Open4::popen4 "R --vanilla --slave"
    stdin.write "source('#{File.join(Sent.rdir,'matrix.R')}');\n"
    stdin.write "#{ command };\n"   
    stdin.close

    Process.wait pid
    raise ProcessAborted, "Error in R process" if $?.exitstatus < 0
    result = stdout.read + stderr.read
    stdout.close
    stderr.close

    result
  end

  def self.metadocs(assocfile, output, low=0.001, hi=0.65, max=3000)
    require 'rbbt/bow/bow'
    require 'rbbt/bow/dictionary'
    
    associations = Open.to_hash(assocfile, :flatten => true, :sep => "\t|,")

    dict = Dictionary::TF_IDF.new

    String.reset_stem_list

    require 'progress-meter'
    Progress.monitor("Building Dictionary for #{File.basename(output)}", 1000) 
    associations.each{|gene, pmids|
      text = PubMed.get_article(pmids).collect{|p| p[1].text}.join("\n")
      dict.add(BagOfWords.count(text.bigrams))
    }

    # At least 3 genes must have a word to be chosen
    hard_min = 3 * 100 / associations.keys.length
    hi = hard_min if hi < hard_min
    
    d =  dict.weights(:low => low, :hi => hi, :limit => max)
    Open.write(output + '.dict', d.sort.collect{|p| p.join("\t")}.join("\n"))
    terms = d.keys.sort

    fout = File.open(output, 'w')
    fout.puts("\t" + terms.join("\t"))

    Progress.monitor("Building Metadoc for #{File.basename(output)}", 1000) 
    associations.each{|gene, pmids|
      text = PubMed.get_article(pmids).collect{|p| p[1].text}.join("\n")
      fout.puts(([gene] + BagOfWords.features(text, terms)).join("\t"))
    }
    fout.close

    Open.write(output + '.stems', String.stem_list(terms.collect{|p| p.split(/ /)}.flatten.uniq).collect{|k,v| "#{ k }\t#{v.join("\t")}"}.join("\n"))
  end

  def self.matrix(metadocs, output, list=nil)
    list ||= []
      
    if list.empty?
      FileUtils.cp metadocs, output
    else
      `head -n 1 #{ metadocs } > #{ output }`
      `grep '^\\(#{list.join('\\|')}\\)[[:space:]]' #{ metadocs } >> #{output}`
      raise Sent::NoGenes, "No Genes Matched" if $? != 0
    end

    dict = metadocs + '.dict'
    run_R("SENT.prepare.matrix('#{ output }', '#{ output }', '#{metadocs + '.dict'}')")
  end


  @@bionmf_wsdl = "http://bionmf.dacya.ucm.es/WebService/BioNMFWS.wsdl"   
  def self.NMF(matrix, out, k, executions = 10)
    driver = SOAP::WSDLDriverFactory.new( @@bionmf_wsdl).create_rpc_driver

    # Upload matrix
    nmf_matrix = driver.upload_matrix(
       File.open(matrix).read, # matrix
       false, # binary 
       true,  # column labels 
       true,  # row labels
       true,  # transpose
       "No",  # positive
       "No",  # normalization 
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
            k, # Factor Start
            k, # Factor End
            1, # Runs
            2000, # Iterations
            40,   # Stop criteria
            0,    # Not used (nsnmf smoothness)
            false, # extra info
            '')    # Suggested name

          while !driver.done(job_id)
            sleep(5)
          end

          if driver.error(job_id) 
            error = driver.messages(job_id).join("\n")
            raise "Error in NMF"  + driver.messages(job_id).join("\n")
          end

          results =  driver.results(job_id)
          fw = File.open(out + ".matrix_w.#{num}",'w')
          fw.write(driver.result(results[0]).sub(/\t(.*)\t$/,'\1'))
          fw.close
          fh = File.open(out + ".matrix_h.#{num}",'w')
          fh.write(driver.result(results[1]).sub(/\t(.*)\t$/,'\1'))
          fh.close
          driver.clean(job_id)
        rescue Sent::ProcessAborted
          puts "Process aborted for #{ num }"
          driver.abort(job_id)
        rescue Timeout::Error
          if times > 0
            times -= 1
            sleep 2
            retry
          else
            raise Sent::ProcessAborted, "NMF Execution #{ num } timed out"
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
            raise Sent::ProcessAborted, "NMF Execution #{ num } Produced Exception"
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
      threads.each{|t| t.raise Sent::ProcessAborted, "Process Aborted"}
      aborted = true
    end

    threads.each { |aThread|  aThread.join }

    Signal.trap("INT", old_int)
    driver.clean(nmf_matrix)

    if aborted
      raise Sent::ProcessAborted, "Process Aborted"
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
      raise Sent::ProcessAborted, "Process Aborted"
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
    run_R("SENT.analyze('#{ prefix }', '#{ output }', '#{clusters}', '#{num_words}')")
    words = Dir.glob(output + '*.words').sort.collect{|f| Open.read(f).split(/\n/)}
    genes = Dir.glob(output + '*.genes').sort.collect{|f| Open.read(f).split(/\n/)}

    groups = []
    words.zip(genes).each{|p|
      groups << {:words => p[0], :genes => p[1]}
    }

    require 'yaml'
    Open.write(output + '.summary', groups.to_yaml)
  end

  def self.literature_index(pmids, outfile)

    index = Ferret::Index::Index.new(:path => outfile)

    index.field_infos.add_field(:title,  :index => :yes, :boost => 0.67)
    index.field_infos.add_field(:abstract, :index => :yes, :boost => 0.33)

    require 'progress-meter'
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