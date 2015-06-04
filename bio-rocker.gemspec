Gem::Specification.new do |s|
   s.name	= 'bio-rocker'
   s.version	= '0.2.13'
   s.license	= 'artistic 2.0'
   s.summary	= 'ROCker'
   s.description = 'Detecting and quantifying functional genes in short-read metagenomic datasets'
   s.authors	= ["Luis (Coto) Orellana","Luis M. Rodriguez-R"]
   s.email	= "lhorellana@gatech.edu"
   s.files	= ["lib/rocker.rb", "lib/rocker/sequence.rb", "lib/rocker/alignment.rb", "lib/rocker/blasthit.rb", "lib/rocker/rocwindow.rb", "lib/rocker/rocdata.rb", "lib/rocker/rinterface.rb"]
   s.homepage	= 'http://enve-omics.ce.gatech.edu/rocker'
   s.executables << 'ROCker'
   s.date	= '2015-06-04'
   s.add_runtime_dependency 'rest-client', '~> 1.7.3'
   s.add_runtime_dependency 'json'
end

