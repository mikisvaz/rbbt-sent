$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require 'yaml'
require 'fileutils'
require 'stemmer'
require 'ferret'

module Sent

  class NoConfigError < StandardError; end

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



