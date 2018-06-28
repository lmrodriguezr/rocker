$:.push File.expand_path("lib", File.dirname(__FILE__))
require "rocker"

Gem::Specification.new do |s|
   s.name	= "bio-rocker"
   s.version	= ROCker.VERSION
   s.license	= "Artistic-2.0"
   s.summary	= "ROCker"
   s.description = "Detecting and quantifying functional genes in " +
		     "short-read metagenomic datasets"
   s.authors	= ["Luis (Coto) Orellana","Luis M. Rodriguez-R"]
   s.email	= "lhorellana@gatech.edu"
   s.files	= ["lib/rocker.rb"]
   s.files	+= %w{sequence.rb alignment.rb blasthit.rb rocwindow.rb
		     rocdata.rb rinterface.rb genome-set.rb
		     protein-set.rb}.map{|f| "lib/rocker/#{f}"}
   s.files	+= %w{build.rb compile.rb search.rb filter.rb
		     plot.rb}.map{|f| "lib/rocker/step/#{f}"}
   s.homepage	= "http://enve-omics.ce.gatech.edu/rocker"
   s.executables << "ROCker"
   s.date	= ROCker.DATE
   s.add_runtime_dependency "rest-client", "~> 1.7"
   s.add_runtime_dependency "json", "~> 1.8"
   s.required_ruby_version = ">= 2.0"
end

