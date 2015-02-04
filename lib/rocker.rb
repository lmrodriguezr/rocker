#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Feb-04-2015
#

require 'rocker/blasthit'
require 'rocker/rocdata'

class ROCker
   #================================[ Class ]
   @@EUTILS = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
   @@DEFAULTS = {
      # General
      :q=>false, :r=>'R', :nucl=>false, :debug=>false,
      # Build
      :positive=>[], :negative=>[], :thr=>2,:genomefrx=>1.0, :pergenus=>false, :perspecies=>false,
	 # ext. software
	 :grinder=>'grinder', :muscle=>'muscle', :blastbins=>'', :seqdepth=>3, :minovl=>0.75,
	 :grindercmd=>'%1$s -reference_file "%2$s" -cf "%3$f" -base_name "%4$s" -dc \'-~*Nn\' -md "poly4 3e-3 3.3e-8" -mr "95 5" -rd "100 uniform 5"',
	 :musclecmd=>'%1$s -in "%2$s" -out "%3$s" -quiet',
	 :blastcmd=>'%1$s%2$s -query "%3$s" -db "%4$s" -out "%5$s" -num_threads %6$d -outfmt 6 -max_target_seqs 1',
	 :makedbcmd=>'%1$smakeblastdb -dbtype %2$s -in "%3$s" -out "%4$s"',
      # Compile
      :refine=>true, :win=>20, :minscore=>0,
      # Filter
      :sbj=>[],
      # Plot
      :color=>false, :gformat=>'pdf', :width=>9, :height=>9, :impact=>false, :transparency=>true,
   }
   @@HAS_BUILD_GEMS = nil
   def self.eutils() @@EUTILS end
   def self.defaults() @@DEFAULTS end
   def self.default(k) @@DEFAULTS[k] end
   def self.has_build_gems?
      return @@HAS_BUILD_GEMS unless @@HAS_BUILD_GEMS.nil?
      @@HAS_BUILD_GEMS = TRUE
      begin
	 require 'rubygems'
	 require 'restclient'
	 require 'nokogiri'
      rescue LoadError
	 @@HAS_BUILD_GEMS = FALSE
      end
      @@HAS_BUILD_GEMS
   end

   #================================[ Instance ]
   attr_reader :o
   def initialize(opts)
      @o = ROCker.defaults
      opts.each{ |k,v| @o[k] = v }
      RInterface.R_BIN = opts[:r] unless opts[:r].nil?
   end
   
   #================================[ Build ]
   def build!
      # Check requirements
      puts "Testing environment." unless @o[:q]
      @o[:noblast]=true if @o[:nomg]
      raise "Unsatisfied requirements, please see the help message (-h)." unless ROCker.has_build_gems?
      @o[:positive] += @o[:posori] unless @o[:posori].nil?
      @o[:positive] += File.readlines(@o[:posfile]).map{ |l| l.chomp } unless @o[:posfile].nil?
      @o[:negative] += File.readlines(@o[:negfile]).map{ |l| l.chomp } unless @o[:negfile].nil?
      unless @o[:aln].nil?
         aln = Alignment.new
	 aln.read_fasta @o[:aln]
	 @o[:positive] += aln.get_gis
      end
      raise "-p or -P are mandatory." if @o[:positive].size==0
      raise "-o/--baseout is mandatory." if @o[:baseout].nil?
      if @o[:positive].size == 1 and not @o[:noaln]
	 warn "\nWARNING: Positive set contains only one sequence, turning off alignment.\n\n"
	 @o[:noaln] = true
      end
      self.bash "#{@o[:grinder]} --version", "-G/--grinder must be executable. Is Grinder installed?" unless @o[:nomg]
      self.bash "#{@o[:muscle]} -version", "-M/--muscle must be executable. Is Muscle installed?" unless @o[:noaln]
      self.bash "#{@o[:blastbins]}makeblastdb -version", "-B/--blastbins must contain executables. Is BLAST+ installed?" unless @o[:noblast]
      # Download genes
      puts "Downloading gene data." unless @o[:q]
      f = File.open(@o[:baseout] + '.ref.fasta', 'w')
      if @o[:posori].nil? and @o[:posfile].nil? and not @o[:aln].nil?
	 puts "  * re-using aligned sequences as positive set." unless @o[:q]
	 f.print aln.to_seq_s
	 @o[:noaln] = true
      else
	 puts "  * downloading #{@o[:positive].size} sequence(s) in positive set." unless @o[:q]
	 $stderr.puts "   # #{@o[:positive]}" if @o[:debug]
	 ids = Array.new(@o[:positive])
	 while ids.size>0
	    f.print efetch({:db=>(@o[:nucl] ? 'nuccore' : 'protein'), :id=>ids.shift(200).join(','), :rettype=>'fasta', :retmode=>'text'})
	 end
      end
      f.close
      genome_gis = {:positive=>[], :negative=>[]}
      [:positive, :negative].each do |set|
         unless @o[set].size==0
	    puts "  * gathering genomes from #{@o[set].size} #{set.to_s} sequence(s)." unless @o[:q]
	    $stderr.puts "   # #{@o[set]}" if @o[:debug]
	    genome_gis[set] = genes2genomes(@o[set], @o[:nucl])
	 end
      end
      raise "No genomes associated with the positive set." if genome_gis[:positive].size==0
      genome_gis[:positive] = genome_gis[:positive].sample( (genome_gis[:positive].size*@o[:genomefrx]).round ) if @o[:genomefrx]
      raise "No positive genomes selected for metagenome construction, is --genome-frx too small?" if genome_gis[:positive].empty?
      all_gis = genome_gis.values.reduce(:+).uniq
      
      # Locate genes
      puts "Analyzing genome data." unless @o[:q]
      puts "  * downloading and parsing #{genome_gis[:positive].size} XML file(s)." unless @o[:q]
      $stderr.puts "   # #{genome_gis[:positive]}" if @o[:debug]
      positive_coords = {}
      genome_org = {}
      i = 0
      genome_gis[:positive].each do |gi|
	 print "  * scanning #{(i+=1).ordinalize} genome out of #{genome_gis[:positive].size}. \r" unless @o[:q]
	 $stderr.puts "   # Looking for any of #{@o[:positive]}" if @o[:debug]
	 genome_file = @o[:baseout] + '.src.' + i.to_s + '.xml'
	 if @o[:reuse] and File.exist? genome_file
	    puts "  * reusing existing file: #{genome_file}." unless @o[:q]
	    ifh = File.open(genome_file, 'r')
	    doc = Nokogiri::XML( ifh )
	    ifh.close
	 else
	    genome_file=nil unless @o[:noclean]
	    res = efetch({:db=>'nuccore', :id=>gi, :rettype=>'xml', :retmode=>'text'}, genome_file)
	    doc = Nokogiri::XML( res )
	 end
	 incomplete = true
	 doc.xpath('//Bioseq-set/Bioseq-set_seq-set/Seq-entry').each do |genome|
	    genome_gi = genome.at_xpath('./Seq-entry_set/Bioseq-set/Bioseq-set_seq-set/Seq-entry/Seq-entry_seq/Bioseq/Bioseq_id/Seq-id/Seq-id_gi')
	    if !genome_gi.nil? and gi==genome_gi.content
	       incomplete = false
	       positive_coords[gi] ||= []
	       $stderr.puts "\n   # got #{gi}, scanning" if @o[:debug]
	       if @o[:pergenus] or @o[:perspecies]
		  name = genome.at_xpath('./Seq-entry_set/Bioseq-set/Bioseq-set_descr/Seq-descr/Seqdesc/Seqdesc_source/BioSource/BioSource_org/Org-ref/Org-ref_orgname/OrgName/OrgName_name/OrgName_name_binomial/BinomialOrgName')
		  unless name.nil?
		     name_g = name.at_xpath('./BinomialOrgName_genus')
		     name_s = name.at_xpath('./BinomialOrgName_species')
		     if name_g.nil? or (name_s.nil? and @o[:perspecies])
		        name = nil
		     else
		        name = @o[:perspecies] ? name_g.content + " " + name_s.content :  name_g.content
		     end
		  end
		  if name.nil?
		     warn "WARNING: Cannot find binomial name of #{gi}, using genome regardless of taxonomy."
		     name = rand(36**100).to_s(36)
		  end
		  break unless genome_org[ name ].nil?
		  genome_org[ name ] = gi
	       end
	       $stderr.puts "   # traversing #{gi}" if @o[:debug]
	       genome.xpath('./Seq-entry_set/Bioseq-set/Bioseq-set_annot/Seq-annot/Seq-annot_data/Seq-annot_data_ftable/Seq-feat').each do |pr|
		  pr_gi = pr.at_xpath('./Seq-feat_product/Seq-loc/Seq-loc_whole/Seq-id/Seq-id_gi')
		  next if pr_gi.nil?
		  if @o[:positive].include? pr_gi.content
		     $stderr.puts "   # found #{pr_gi.content}" if @o[:debug]
		     pr_loc = pr.at_xpath('./Seq-feat_location/Seq-loc/Seq-loc_int/Seq-interval')
		     if pr_loc.nil?
			pr_loc = pr.xpath('./Seq-feat_location/Seq-loc/Seq-loc_mix//Seq-loc/Seq-loc_int/Seq-interval')
			if pr_loc.nil?
			   warn "WARNING: Impossible to find location of '#{pr_gi.content}' in '#{gi}'."
			   incomplete = true
			else
			   pr_loc.each do |loc_int|
			      positive_coords[gi] << {
				 :gi     => pr_gi.content,
				 :from   => loc_int.at_xpath('./Seq-interval_from').content.to_i,
				 :to     => loc_int.at_xpath('./Seq-interval_to').content.to_i
				 #, :strand => loc_int.at_xpath('./Seq-interval_strand/Na-strand/@value').content
			      }
			   end
			end
		     else
			positive_coords[gi] << {
			   :gi     => pr_gi.content,
			   :from   => pr_loc.at_xpath('./Seq-interval_from').content.to_i,
			   :to     => pr_loc.at_xpath('./Seq-interval_to').content.to_i
			   #, :strand => pr_loc.at_xpath('./Seq-interval_strand/Na-strand/@value').content
			}
		     end
		  end
	       end
	       break
	    end
	 end
	 doc = nil
	 warn "WARNING: Cannot find GI '#{gi}'." if incomplete
      end
      genome_gis[:positive] = genome_org.values if @o[:pergenus] or @o[:perspecies]
      all_gis = genome_gis.values.reduce(:+).uniq
      print "\n" unless @o[:q]
      missing = @o[:positive] - positive_coords.values.map{ |a| a.map{ |b| b[:gi] } }.reduce(:+)
      warn "\nWARNING: Cannot find genomic location of sequence(s) #{missing.join(',')}.\n\n" unless missing.size==0 or @o[:genomefrx]<1.0 or @o[:pergenus] or @o[:perspecies]
      
      # Download genomes
      genomes_file = @o[:baseout] + '.src.fasta'
      if @o[:reuse] and File.exist? genomes_file
	 puts "  * reusing existing file: #{genomes_file}." unless @o[:q]
      else
	 puts "  * downloading #{all_gis.size} genome(s) in FastA." unless @o[:q]
	 $stderr.puts "   # #{all_gis}" if @o[:debug]
	 ids = Array.new(all_gis)
	 ofh = File.open(genomes_file, 'w')
	 while ids.size>0
	    ofh.print efetch({:db=>'nuccore', :id=>ids.shift(200).join(','), :rettype=>'fasta', :retmode=>'text'})
	 end
	 ofh.close
      end
      
      # Generate metagenome
      unless @o[:nomg]
	 puts "Generating in silico metagenome" unless @o[:q]
	 if @o[:reuse] and File.exist? @o[:baseout] + ".mg.fasta"
	    puts "  * reusing existing file: #{@o[:baseout]}.mg.fasta." unless @o[:q]
	 else
	    all_src = File.readlines("#{@o[:baseout]}.src.fasta").select{ |l| l =~ /^>/ }.size
	    thrs = [@o[:thr], all_src].min
	    puts "  * running grinder and tagging positive reads (#{thrs} threads)." unless @o[:q]
	    $stderr.puts "   # #{positive_coords}" if @o[:debug]
	    thr_obj = []
	    seqs_per_thr = (all_src/thrs).ceil
	    (0 .. (thrs-1)).each do |thr_i|
	       thr_obj << Thread.new do
		  Thread.current[:seqs_a] = thr_i*seqs_per_thr + 1
		  Thread.current[:seqs_b] = [Thread.current[:seqs_a] + seqs_per_thr, all_src].min
		  # Create sub-fasta
		  Thread.current[:ofh] = File.open("#{@o[:baseout]}.src.fasta.#{thr_i.to_s}", 'w')
		  Thread.current[:ifh] = File.open("#{@o[:baseout]}.src.fasta", 'r')
		  Thread.current[:seq_i] = 0
		  while Thread.current[:l] = Thread.current[:ifh].gets
		     Thread.current[:seq_i]+=1 if Thread.current[:l] =~ /^>/
		     break if Thread.current[:seq_i] > Thread.current[:seqs_b]
		     Thread.current[:ofh].print Thread.current[:l] if Thread.current[:seq_i] >= Thread.current[:seqs_a]
		  end
		  Thread.current[:ifh].close
		  Thread.current[:ofh].close
		  bash sprintf(@o[:grindercmd], @o[:grinder], "#{@o[:baseout]}.src.fasta.#{thr_i.to_s}", @o[:seqdepth], "#{@o[:baseout]}.mg.tmp.#{thr_i.to_s}")
		  # Tag positives
		  puts "  * tagging positive reads." unless @o[:q]
		  Thread.current[:ifh] = File.open(@o[:baseout] + ".mg.tmp.#{thr_i.to_s}-reads.fa", 'r')
		  Thread.current[:ofh] = File.open(@o[:baseout] + ".mg.fasta.#{thr_i.to_s}", 'w')
		  while Thread.current[:l]=Thread.current[:ifh].gets
		     Thread.current[:rd] = /^>(?<id>\d+) reference=gi\|(?<gi>\d+)\|.* position=(?<comp>complement\()?(?<from>\d+)\.\.(?<to>\d+)\)? /.match(Thread.current[:l])
		     unless Thread.current[:rd].nil?
			Thread.current[:positive] = false
			positive_coords[Thread.current[:rd][:gi]] ||= []
			positive_coords[Thread.current[:rd][:gi]].each do |gn|
			   Thread.current[:left]  = Thread.current[:rd][:to].to_i - gn[:from]
			   Thread.current[:right] = gn[:to] - Thread.current[:rd][:from].to_i
			   if (Thread.current[:left]*Thread.current[:right] >= 0) and ([Thread.current[:left], Thread.current[:right]].min/(Thread.current[:rd][:to].to_i-Thread.current[:rd][:from].to_i) >= @o[:minovl])
			      Thread.current[:positive] = true
			      break
			   end
			end
			Thread.current[:l] = ">#{Thread.current[:rd][:id]}#{Thread.current[:positive] ? "@%" : ""} ref=#{Thread.current[:rd][:gi]}:#{Thread.current[:rd][:from]}..#{Thread.current[:rd][:to]}#{(Thread.current[:rd][:comp]=='complement(')?'-':'+'}\n"
		     end
		     Thread.current[:ofh].print Thread.current[:l]
		  end
		  Thread.current[:ofh].close
		  Thread.current[:ifh].close
		  Thread.current[:output] = @o[:baseout] + ".mg.fasta.#{thr_i.to_s}"
	       end # Thread.new do
	    end # (1 .. thrs).each
	    # Concatenate results
	    ofh = File.open(@o[:baseout] + ".mg.fasta", 'w')
	    thr_obj.each do |t|
	       t.join
	       raise "Thread failed without error trace: #{t}" if t[:output].nil?
	       ifh = File.open(t[:output], 'r')
	       while l = ifh.gets
	          ofh.print l
	       end
	       ifh.close
	       File.unlink t[:output]
	    end
	    ofh.close
         end
      end # unless @o[:nomg]
      # Align references
      unless @o[:noaln]
	 puts "Aligning reference set." unless @o[:q]
	 if @o[:reuse] and File.exist? "#{@o[:baseout]}.ref.aln"
	    puts "  * reusing existing file: #{@o[:baseout]}.ref.aln." unless @o[:q]
	 else
	    bash sprintf(@o[:musclecmd], @o[:muscle], "#{@o[:baseout]}.ref.fasta", "#{@o[:baseout]}.ref.aln")
	    puts "  +--\n  | IMPORTANT NOTE: Manually checking the alignment before\n  | the 'compile' step is *strongly* encouraged.\n  +--\n" unless @o[:q]
	 end
      end
      # Run BLAST 
      unless @o[:noblast]
	 puts "Running homology search." unless @o[:q]
	 if @o[:reuse] and File.exist? "#{@o[:baseout]}.ref.blast"
	    puts "  * reusing existing file: #{@o[:baseout]}.ref.blast." unless @o[:q]
	 else
	    puts "  * preparing database." unless @o[:q]
	    bash sprintf(@o[:makedbcmd], @o[:blastbins], (@o[:nucl]?'nucl':'prot'), "#{@o[:baseout]}.ref.fasta", "#{@o[:baseout]}.ref")
	    puts "  * running BLAST." unless @o[:q]
	    bash sprintf(@o[:blastcmd], @o[:blastbins], (@o[:nucl]?'blastn':'blastx'), "#{@o[:baseout]}.mg.fasta", "#{@o[:baseout]}.ref", "#{@o[:baseout]}.ref.blast", @o[:thr])
	 end
      end
      # Clean
      unless @o[:noclean]
	 puts "Cleaning." unless @o[:q]
	 sff  = %w{.src.xml .src.fasta}
	 sff += %w{.mg.tmp-reads.fa .mg.tmp-ranks.txt} unless @o[:nomg]
	 sff += %w{.ref.phr .ref.pin .ref.psq} unless @o[:noblast]
	 sff.each { |sf| File.unlink @o[:baseout] + sf if File.exist? @o[:baseout] + sf }
      end
   end # build!

   #================================[ Compile ]
   def compile!
      raise "-a/--alignment is mandatory." if @o[:aln].nil?
      raise "-a/--alignment must exist." unless File.exist? @o[:aln]
      if @o[:table].nil?
	 raise "-t/--table is mandatory unless -b is provided." if @o[:blast].nil?
	 @o[:table] = "#{@o[:blast]}.table"
      end
      raise "-b/--blast is mandatory unless -t exists." if @o[:blast].nil? and not File.exist? @o[:table]
      raise "-k/--rocker is mandatory." if @o[:rocker].nil?

      puts "Testing environment." unless @o[:q]
      bash "echo '' | #{@o[:r]} --vanilla", "-r/--path-to-r must be executable. Is R installed?"
      bash "echo \"library('pROC')\" | #{@o[:r]} --vanilla", "Please install the 'pROC' library for R first."

      puts "Reading files." unless @o[:q]
      puts "  * loading alignment: #{@o[:aln]}." unless @o[:q]
      aln = Alignment.new
      aln.read_fasta @o[:aln]
      
      if File.exist? @o[:table]
	 puts "  * reusing existing file: #{@o[:table]}." unless @o[:q]
      else
	 puts "  * generating table: #{@o[:table]}." unless @o[:q]
	 blast2table(@o[:blast], @o[:table], aln, @o[:minscore])
      end

      puts "Analyzing data." unless @o[:q]
      puts "  * computing windows." unless @o[:q]
      data = ROCData.new(@o[:table], aln, @o[:win])
      data.nucl = @o[:nucl]
      if @o[:refine]
	 puts "  * refining windows." unless @o[:q]
	 warn "Insufficient hits to refine results." unless data.refine! @o[:table]
      end
      puts "  * saving ROCker file: #{@o[:rocker]}." unless @o[:q]
      data.save @o[:rocker]
   end # compile!
   
   #================================[ Filter ]
   def filter!
      raise "-k/--rocker is mandatory." if @o[:rocker].nil?
      raise "-x/--query-blast is mandatory." if @o[:qblast].nil?
      raise "-o/--out-blast is mandatory." if @o[:oblast].nil?
      
      puts "Reading ROCker file." unless @o[:q]
      data = ROCData.new @o[:rocker]

      puts "Filtering BLAST." unless @o[:q]
      ih = File.open(@o[:qblast], 'r')
      oh = File.open(@o[:oblast], 'w')
      while ln = ih.gets
	 bh = BlastHit.new(ln, data.aln)
	 oh.print ln if not(bh.sfrom.nil?) and bh.bits >= data.win_at_col(bh.midpoint).thr
      end
      ih.close
      oh.close
   end # filter!
   #================================[ Search ]
   def search!
      raise "-k/--rocker is mandatory." if @o[:rocker].nil?
      raise "Code Under development..."
      # ToDo
      # [ ... ]
   end # search!
   
   #================================[ Plot ]
   def plot!
      raise "-k/--rocker is mandatory." if o[:rocker].nil?
      if @o[:table].nil?
	 raise "-t/--table is mandatory unless -b is provided." if @o[:blast].nil?
	 @o[:table] = "#{@o[:blast]}.table"
      end
      raise "-b/--blast is mandatory unless -t exists." if @o[:blast].nil? and not File.exist? @o[:table]

      puts "Testing environment." unless @o[:q]
      bash "echo '' | #{@o[:r]} --vanilla", "-r/--path-to-r must be executable. Is R installed?"

      puts "Reading files." unless @o[:q]
      puts "  * loding ROCker file: #{@o[:rocker]}." unless @o[:q]
      data = ROCData.new @o[:rocker]
      if File.exist? @o[:table]
	 puts "  * reusing existing file: #{@o[:table]}." unless @o[:q]
      else
	 puts "  * generating table: #{@o[:table]}." unless @o[:q]
	 blast2table(@o[:blast], @o[:table], data.aln, @o[:minscore])
      end

      puts "Plotting hits." unless @o[:q]
      extra = @o[:gformat]=='pdf' ? "" : ", units='in', res=300"
      @o[:gout] ||= "#{@o[:rocker]}.#{@o[:gformat]}"
      data.rrun "#{@o[:gformat]}('#{@o[:gout]}', #{@o[:width]}, #{@o[:height]}#{extra});"
      data.rrun "layout(c(2,1,3), heights=c(2-1/#{data.aln.size},3,1));"
      some_thr = data.load_table! @o[:table], @o[:sbj], @o[:minscore]
      data.rrun "par(mar=c(0,4,0,0.5)+.1);"
      data.rrun "plot(1, t='n', xlim=c(0.5,#{data.aln.cols}+0.5), ylim=range(x$V4)+c(-0.04,0.04)*diff(range(x$V4)), xlab='', ylab='Bit score', xaxs='i', xaxt='n');"
      data.rrun "noise <- runif(ncol(x),-.2,.2)"
      data.rrun "arrows(x0=x$V2, x1=x$V3, y0=x$V4+noise, col=ifelse(x$V5==1, rgb(0,0,.5,#{@o[:transparency] ? ".2" : "1"}), rgb(.5,0,0,#{@o[:transparency] ? ".2" : "1"})), length=0);"
      data.rrun "points(x$V6, x$V4+noise, col=ifelse(x$V5==1, rgb(0,0,.5,#{@o[:transparency] ? ".5" : "1"}), rgb(.5,0,0,#{@o[:transparency] ? ".5" : "1"})), pch=19, cex=1/4);"

      puts "Plotting windows." unless @o[:q]
      if some_thr
	 data.rrun "arrows(x0=w$V1, x1=w$V2, y0=w$V5, lwd=2, length=0)"
	 data.rrun "arrows(x0=w$V2[-nrow(w)], x1=w$V1[-1], y0=w$V5[-nrow(w)], y1=w$V5[-1], lwd=2, length=0)"
      end
      data.rrun "legend('bottomright',legend=c('Hit span','Hit mid-point','Reference','Non-reference')," +
	 "lwd=c(1,NA,1,1),pch=c(NA,19,19,19),col=c('black','black','darkblue','darkred'),ncol=4,bty='n')"

      puts "Plotting alignment." unless @o[:q]
      data.rrun "par(mar=c(0,4,0.5,0.5)+0.1);"
      data.rrun "plot(1, t='n', xlim=c(0,#{data.aln.cols}),ylim=c(1,#{data.aln.seqs.size}),xlab='',ylab='Alignment',xaxs='i',xaxt='n',yaxs='i',yaxt='n',bty='n');"
      i = 0
      data.rrun "clr <- rainbow(26, v=1/2, s=3/4);" if @o[:color]
      data.aln.seqs.values.each do |s|
         color = s.aln.split(//).map{|c| c=="-" ? "'grey80'" : (@o[:sbj].include?(s.id) ? "'red'" : (@o[:color] ? "clr[#{c.ord-64}]" : "'black'"))}.join(',')
	 data.rrun "rect((1:#{data.aln.cols-1})-0.5, rep(#{i}, #{data.aln.cols-1}), (1:#{data.aln.cols-1})+0.5, rep(#{i+1}, #{data.aln.cols-1}), col=c(#{color}), border=NA);"
	 i += 1
      end

      puts "Plotting statistics." unless @o[:q]
      data.rrun "par(mar=c(5,4,0,0.5)+.1);"
      data.rrun "plot(1, t='n', xlim=c(0,#{data.aln.cols}),ylim=c(#{@o[:impact] ? "-2,.1" : "50,100"}),xlab='Alignment position (amino acids)',ylab='Precision',xaxs='i');"
      if some_thr
	 sn = data.rrun "100*sum(w$tp)/(sum(w$tp)+sum(w$fn))", :float
	 sp = data.rrun "100*sum(w$tn)/(sum(w$fp)+sum(w$tn))", :float
	 ac = data.rrun "100*(sum(w$tp)+sum(w$tn))/(sum(w$p)+sum(w$n))", :float
	 unless @o[:q]
	    puts "  * sensitivity: #{sn}%"
	    puts "  * specificity: #{sp}%"
	    puts "  * accuracy: #{ac}%"
	 end
	 data.rrun "pos <- (w$V1+w$V2)/2"
	 if @o[:impact]
	    data.rrun "lines(pos[!is.na(w$specificity)], (w$specificity[!is.na(w$specificity)]-#{sp})*w$tp[!is.na(w$specificity)]/sum(w$tp), col='darkred', lwd=2, t='o', cex=1/3, pch=19);"
	    data.rrun "lines(pos[!is.na(w$sensitivity)], (w$sensitivity[!is.na(w$sensitivity)]-#{sn})*w$tn[!is.na(w$sensitivity)]/sum(w$tn), col='darkgreen', lwd=2, t='o', cex=1/3, pch=19);"
	    data.rrun "lines(pos[!is.na(w$accuracy)], (w$accuracy[!is.na(w$accuracy)]-#{ac})*(w$tp+w$tn)[!is.na(w$accuracy)]/sum(c(w$tp, w$tn)), col='darkblue', lwd=2, t='o', cex=1/3, pch=19);"
	 else
	    data.rrun "lines(pos[!is.na(w$specificity)], w$specificity[!is.na(w$specificity)], col='darkred', lwd=2, t='o', cex=1/3, pch=19);"
	    data.rrun "lines(pos[!is.na(w$sensitivity)], w$sensitivity[!is.na(w$sensitivity)], col='darkgreen', lwd=2, t='o', cex=1/3, pch=19);"
	    data.rrun "lines(pos[!is.na(w$accuracy)], w$accuracy[!is.na(w$accuracy)], col='darkblue', lwd=2, t='o', cex=1/3, pch=19);"
	 end
	 #data.rrun "lines(pos[!is.na(w$precision)], w$precision[!is.na(w$precision)], col='purple', lwd=2, t='o', cex=1/3, pch=19);"
      end
      data.rrun "legend('bottomright',legend=c('Specificity','Sensitivity','Accuracy'),lwd=2,col=c('darkred','darkgreen','darkblue'),ncol=3,bty='n')"
      data.rrun "dev.off();"
   end # plot!

   #================================[ Utilities ]
   def blast2table(blast_f, table_f, aln, minscore)
      ifh = File.open(blast_f, "r")
      ofh = File.open(table_f, "w")
      while ln = ifh.gets
	 bh = BlastHit.new(ln, aln)
	 ofh.print bh.to_s if bh.bits >= minscore
      end
      ifh.close
      ofh.close
   end
   def genes2genomes(gis, nucl=false)
      genomes = []
      ids = Array.new(gis)
      while ids.size>0
	 doc = Nokogiri::XML( elink({:dbfrom=>(nucl ? 'nuccore' : 'protein'), :db=>'nuccore', :id=>ids.shift(200).join(',')}) )
	 genomes += doc.xpath('/eLinkResult/LinkSet/LinkSetDb/Link/Id').map{ |id| id.content }
      end
      genomes.uniq
   end
   def eutils(script, params={}, outfile=nil)
      response = RestClient.get "#{ROCker.eutils}/#{script}", {:params=>params}
      raise "Unable to reach NCBI EUtils, error code #{response.code}." unless response.code == 200
      unless outfile.nil?
	 ohf = File.open(outfile, 'w')
	 ohf.print response.to_s
	 ohf.close
      end
      response.to_s
   end
   def efetch(*etc) self.eutils 'efetch.fcgi', *etc end
   def elink(*etc) self.eutils 'elink.fcgi', *etc end
   def bash(cmd, err_msg=nil)
      o = `#{cmd} 2>&1 && echo '{'`
      raise (err_msg.nil? ? "Error executing: #{cmd}\n\n#{o}" : err_msg) unless o[-2]=='{'
      true
   end
end

#================================[ Extensions ]
class Numeric
   def ordinalize
      n= self.to_s
      s= n[-2]=='1' ? 'th' :
	 n[-1]=='1' ? 'st' :
	 n[-1]=='2' ? 'nd' :
	 n[-1]=='3' ? 'rd' : 'th'
      n + s
   end
end

