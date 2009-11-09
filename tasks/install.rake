require 'sent'

$datadir = Sent.datadir
$scriptdir = File.join(Sent.rootdir, '/install_scripts')

task 'analysis' do
  directory = "#{$datadir}/analysis"
  FileUtils.mkdir_p directory
  %w(Rakefile).each{|f|
    FileUtils.cp_r File.join($scriptdir, "analysis/#{ f }"), directory 
  }
  
  %w(associations metadocs matrices NMF summary).each{|d|
    FileUtils.mkdir_p File.join(directory, d)
  }
end

