Gem::Specification.new do |s|
   s.name	= 'bio-rocker'
   s.version	= '1.0.6'
   s.license	= 'artistic 2.0'
   s.summary	= 'ROCker'
   s.description = 'Detecting and quantifying functional genes in short-read metagenomic datasets'
   s.authors	= ["Luis (Coto) Orellana","Luis M. Rodriguez-R"]
   s.email	= "lhorellana@gatech.edu"
   s.files	= ["lib/rocker.rb"]
   s.files	+= %w{sequence.rb alignment.rb blasthit.rb rocwindow.rb rocdata.rb rinterface.rb}.map{|f| "lib/rocker/#{f}"}
   s.files	+= %w{build.rb compile.rb search.rb filter.rb plot.rb}.map{|f| "lib/rocker/step/#{f}"}
   s.homepage	= 'http://enve-omics.ce.gatech.edu/rocker'
   s.executables << 'ROCker'
   s.date	= '2015-06-14'
   s.add_runtime_dependency 'rest-client', '~> 1.7.3'
   s.add_runtime_dependency 'json', '~> 1.8.1'
   s.required_ruby_version = '~> 2.0'
end

