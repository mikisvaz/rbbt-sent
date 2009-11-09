# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run the gemspec command
# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{rbbt-sent}
  s.version = "0.0.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Miguel Vazquez"]
  s.date = %q{2009-11-04}
  s.default_executable = %q{sent_config}
  s.description = %q{Use literature mining to find semantic features to describe clusers of genes}
  s.email = %q{miguel.vazquez@fdi.ucm.es}
  s.executables = ["sent_config"]
  s.extra_rdoc_files = [
    "LICENSE",
     "README.rdoc"
  ]
  s.files = [
    "bin/sent_config",
     "install_scripts/analysis/Rakefile",
     "lib/sent.rb",
     "lib/sent/main.rb",
     "tasks/install.rake",
     "test/helper.rb"
  ]
  s.homepage = %q{http://github.com/mikisvaz/rbbt-sent}
  s.rdoc_options = ["--charset=UTF-8"]
  s.require_paths = ["lib"]
  s.rubygems_version = %q{1.3.5}
  s.summary = %q{Semantic Features in Text}
  s.test_files = [
    "test/helper.rb"
  ]

  if s.respond_to? :specification_version then
    current_version = Gem::Specification::CURRENT_SPECIFICATION_VERSION
    s.specification_version = 3

    if Gem::Version.new(Gem::RubyGemsVersion) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<rake>, [">= 0"])
      s.add_runtime_dependency(%q<simpleconsole>, [">= 0"])
      s.add_runtime_dependency(%q<rbbt>, [">= 0"])
      s.add_runtime_dependency(%q<ferret>, [">= 0"])
      s.add_runtime_dependency(%q<stemmer>, [">= 0"])
      s.add_runtime_dependency(%q<progress-monitor>, [">= 0"])
      s.add_runtime_dependency(%q<open4>, [">= 0"])
    else
      s.add_dependency(%q<rake>, [">= 0"])
      s.add_dependency(%q<simpleconsole>, [">= 0"])
      s.add_dependency(%q<rbbt>, [">= 0"])
      s.add_dependency(%q<ferret>, [">= 0"])
      s.add_dependency(%q<stemmer>, [">= 0"])
      s.add_dependency(%q<progress-monitor>, [">= 0"])
      s.add_dependency(%q<open4>, [">= 0"])
    end
  else
    s.add_dependency(%q<rake>, [">= 0"])
    s.add_dependency(%q<simpleconsole>, [">= 0"])
    s.add_dependency(%q<rbbt>, [">= 0"])
    s.add_dependency(%q<ferret>, [">= 0"])
    s.add_dependency(%q<stemmer>, [">= 0"])
    s.add_dependency(%q<progress-monitor>, [">= 0"])
    s.add_dependency(%q<open4>, [">= 0"])
  end
end

