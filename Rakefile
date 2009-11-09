require 'rubygems'
require 'rake'

begin
  require 'jeweler'
  Jeweler::Tasks.new do |gem|
    gem.name = "rbbt-sent"
    gem.summary = %Q{Semantic Features in Text}
    gem.description = %Q{Use literature mining to find semantic features to describe clusers of genes}
    gem.email = "miguel.vazquez@fdi.ucm.es"
    gem.homepage = "http://github.com/mikisvaz/rbbt-sent"
    gem.authors = ["Miguel Vazquez"]
  
    gem.files = Dir['lib/**/*', 'test/*', 'bin/sent_config', 'install_scripts/*/*', 'tasks/*', 'R/*']

    gem.add_dependency('rake')
    gem.add_dependency('simpleconsole')
    gem.add_dependency('rbbt')

    gem.add_dependency('ferret')
    gem.add_dependency('stemmer')
    gem.add_dependency('progress-monitor')
    gem.add_dependency('open4')


    # gem is a Gem::Specification... see http://www.rubygems.org/read/chapter/20 for additional settings
  end
rescue LoadError
  puts "Jeweler (or a dependency) not available. Install it with: sudo gem install jeweler"
end

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  test.libs << 'lib' << 'test'
  test.pattern = 'test/**/test_*.rb'
  test.verbose = true
end

begin
  require 'rcov/rcovtask'
  Rcov::RcovTask.new do |test|
    test.libs << 'test'
    test.pattern = 'test/**/test_*.rb'
    test.verbose = true
  end
rescue LoadError
  task :rcov do
    abort "RCov is not available. In order to run rcov, you must: sudo gem install spicycode-rcov"
  end
end

task :test => :check_dependencies

task :default => :test

require 'rake/rdoctask'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "rbbt-sent #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
