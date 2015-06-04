#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jun-04-2015
#

require 'json'
require 'rocker/blasthit'
require 'rocker/rocdata'

class ROCker
   #================================[ Class ]
   @@EBIREST = 'http://www.ebi.ac.uk/Tools'
   @@DEFAULTS = {
      # General
      :q=>false, :r=>'R', :nucl=>false, :debug=>false,
      # Build
      :positive=>[], :negative=>[], :thr=>2,:genomefrx=>1.0,
	 # ext. software
	 :grinder=>'grinder', :muscle=>'muscle', :blastbins=>'', :seqdepth=>0.03, :readlen=>100, :minovl=>50,
	 :grindercmd=>'%1$s -reference_file "%2$s" -cf "%3$f" -dc \'-~*NnKkMmRrYySsWwBbVvHhDdXx\' -md uniform 0.1 -mr 95 5 -rd %4$d uniform 5 -base_name "%5$s"',
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
   def self.ebirest() @@EBIREST ; end
   def self.defaults() @@DEFAULTS ; end
   def self.default(k) @@DEFAULTS[k] ; end
   def self.has_build_gems?
      return @@HAS_BUILD_GEMS unless @@HAS_BUILD_GEMS.nil?
      @@HAS_BUILD_GEMS = TRUE
      begin
	 require 'rubygems'
	 require 'restclient'
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
	 @o[:positive] += aln.get_ids
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
      ref_file = @o[:baseout] + ".ref.fasta"
      if @o[:posori].nil? and @o[:posfile].nil? and not @o[:aln].nil?
	 puts "  * reusing aligned sequences as positive set." unless @o[:q]
	 f = File.open(ref_file, "w")
	 f.print aln.to_seq_s
	 f.close
	 @o[:noaln] = true
      elsif @o[:reuse] and File.size? ref_file
	 puts "  * reusing positive set: #{ref_file}." unless @o[:q]
      else
	 puts "  * downloading #{@o[:positive].size} sequence(s) in positive set." unless @o[:q]
	 $stderr.puts "   # #{@o[:positive]}" if @o[:debug]
	 ids = Array.new(@o[:positive])
	 f = File.open(ref_file, "w")
	 while ids.size>0
	    f.print ebiFetch(:uniprotkb, ids.shift(200), :fasta)
	 end
	 f.close
      end
      genome_ids = {:positive=>[], :negative=>[]}
      [:positive, :negative].each do |set|
         unless @o[set].size==0
	    puts "  * gathering genomes from #{@o[set].size} #{set.to_s} sequence(s)." unless @o[:q]
	    $stderr.puts "   # #{@o[set]}" if @o[:debug]
	    genome_ids[set] = genes2genomes(@o[set])
	 end
      end
      raise "No genomes associated with the positive set." if genome_ids[:positive].size==0
      genome_ids[:positive] = genome_ids[:positive].sample( (genome_ids[:positive].size*@o[:genomefrx]).round ) if @o[:genomefrx]
      raise "No positive genomes selected for metagenome construction, is --genome-frx too small?" if genome_ids[:positive].empty?
      all_genome_ids = genome_ids.values.reduce(:+).uniq
      
      # Locate genes
      puts "Analyzing genome data." unless @o[:q]
      coords_file = @o[:baseout] + ".src.coords"
      if @o[:reuse] and File.size? coords_file
	 puts "  * reusing coordinates: #{coords_file}." unless @o[:q]
	 c = JSON.parse File.read(coords_file), {:symbolize_names=>true}
	 positive_coords = c[:positive_coords]
	 genome_org = c[:genome_org]
      else
	 puts "  * downloading and parsing #{genome_ids[:positive].size} GFF3 document(s)." unless @o[:q]
	 $stderr.puts "   # #{genome_ids[:positive]}" if @o[:debug]
	 positive_coords = {}
	 genome_org = {}
	 i = 0
	 $stderr.puts "   # Looking for any of #{@o[:positive]}" if @o[:debug]
	 genome_ids[:positive].each do |genome_id|
	    print "  * scanning #{(i+=1).ordinalize} genome out of #{genome_ids[:positive].size}. \r" unless @o[:q]
	    unless @o[:pertaxon].nil?
	       genome_taxon = genome2taxon(genome_id, @o[:pertaxon])
	       next unless genome_org[ genome_taxon ].nil?
	       genome_org[ genome_taxon ] = genome_id
	    end
	    genome_file = @o[:baseout] + '.src.' + i.to_s + '.gff3'
	    if @o[:reuse] and File.size? genome_file
	       ifh = File.open(genome_file, 'r')
	       doc = ifh.readlines.grep(/^[^#]/)
	       ifh.close
	    else
	       genome_file=nil unless @o[:noclean]
	       doc = ebiFetch(:embl, [genome_id], :gff3, genome_file).split("\n").grep(/^[^#]/)
	    end
	    doc.each do |ln|
	       next if ln =~ /^#/
	       r = ln.chomp.split /\t/
	       next if r.size < 9
	       prots = r[8].split(/;/).grep(/^db_xref=UniProtKB[\/A-Za-z-]*:/){ |xref| xref.split(/:/)[1] }
	       p = prots.select{ |p| @o[:positive].include? p }.first
	       next if p.nil?
	       positive_coords[ r[0] ] ||= []
	       positive_coords[ r[0] ] << {
		  #:strand	=> r[6],
		  :prot_id	=> p,
		  :from	=> r[3].to_i,
		  :to	=> r[4].to_i
	       }
	    end
	    ofh = File.open(coords_file, 'w')
	    ofh.print JSON.pretty_generate({:positive_coords=>positive_coords, :genome_org=>genome_org})
	    ofh.close
	 end
	 print "\n" unless @o[:q]
      end
      unless @o[:pertaxon].nil?
	 genome_ids[:positive] = genome_org.values
	 puts "  Using #{genome_org.size} genome(s) after filtering by #{@o[:pertaxon]}." unless @o[:q]
      end
      all_genome_ids = genome_ids.values.reduce(:+).uniq
      found = positive_coords.values.map{ |a| a.map{ |b| b[:prot_id] } }.reduce(:+)
      raise "Cannot find the genomic location of any provided sequence." if found.nil?
      missing = @o[:positive] - found
      warn "\nWARNING: Cannot find genomic location of sequence(s) #{missing.join(',')}.\n\n" unless missing.size==0 or @o[:genomefrx]<1.0 or not @o[:pertaxon].nil?
      
      # Download genomes
      genomes_file = @o[:baseout] + '.src.fasta'
      if @o[:reuse] and File.exist? genomes_file
	 puts "  * reusing existing file: #{genomes_file}." unless @o[:q]
      else
	 puts "  * downloading #{all_genome_ids.size} genome(s) in FastA." unless @o[:q]
	 $stderr.puts "   # #{all_genome_ids}" if @o[:debug]
	 ids = Array.new(all_genome_ids)
	 ofh = File.open(genomes_file, 'w')
	 while ids.size>0
	    ofh.print ebiFetch('embl', ids.shift(200), 'fasta')
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

		  # Run grinder (except if the temporal file is already there and can be reused)
		  unless @o[:reuse] and File.size? @o[:baseout] + ".mg.tmp.#{thr_i.to_s}-reads.fa"
		     bash sprintf(@o[:grindercmd], @o[:grinder], "#{@o[:baseout]}.src.fasta.#{thr_i.to_s}", @o[:seqdepth]*@o[:readlen].to_f, @o[:readlen], "#{@o[:baseout]}.mg.tmp.#{thr_i.to_s}")
		  end

		  # Tag positives
		  puts "  * tagging positive reads [thread #{thr_i.to_s}]." unless @o[:q]
		  Thread.current[:ifh] = File.open(@o[:baseout] + ".mg.tmp.#{thr_i.to_s}-reads.fa", 'r')
		  Thread.current[:ofh] = File.open(@o[:baseout] + ".mg.fasta.#{thr_i.to_s}", 'w')
		  while Thread.current[:l]=Thread.current[:ifh].gets
		     if Thread.current[:l] =~ /^>/
			Thread.current[:rd] = /^>(?<id>\d+) reference=[A-Za-z]+\|(?<genome_id>[A-Za-z0-9_]+)\|.* position=(?<comp>complement\()?(?<from>\d+)\.\.(?<to>\d+)\)? /.match(Thread.current[:l])
			raise "Cannot parse simulated read's defline, are you using grinder?: #{Thread.current[:l]}" if Thread.current[:rd].nil?
			Thread.current[:positive] = false
			positive_coords[Thread.current[:rd][:genome_id]] ||= []
			positive_coords[Thread.current[:rd][:genome_id]].each do |gn|
			   Thread.current[:left]  = Thread.current[:rd][:to].to_i - gn[:from]
			   Thread.current[:right] = gn[:to] - Thread.current[:rd][:from].to_i
			   if (Thread.current[:left]*Thread.current[:right] >= 0) and ([Thread.current[:left], Thread.current[:right]].min >= @o[:minovl])
			      Thread.current[:positive] = true
			      break
			   end
			end
			Thread.current[:l] = ">#{thr_i.to_s}_#{Thread.current[:rd][:id]}#{Thread.current[:positive] ? "@%" : ""} " +
			   "ref=#{Thread.current[:rd][:genome_id]}:#{Thread.current[:rd][:from]}..#{Thread.current[:rd][:to]}#{(Thread.current[:rd][:comp]=='complement(')?'-':'+'}\n"
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

      puts "Plotting matches." unless @o[:q]
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
      data.rrun "legend('bottomright',legend=c('Match span','Match mid-point','Reference','Non-reference')," +
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
      data.rrun "plot(1, t='n', xlim=c(0,#{data.aln.cols}),ylim=c(#{@o[:ylim].nil? ? (@o[:impact] ? "-2,.1" : "50,100") : @o[:ylim]}),xlab='Alignment position (amino acids)',ylab='Precision',xaxs='i');"
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
   def genes2genomes(gene_ids)
      genomes = []
      ids = Array.new(gene_ids)
      while ids.size>0
	 doc = ebiFetch(:uniprotkb, ids.shift(200), :annot).split("\n")
	 genomes += doc.grep( /^DR\s+EMBL;/ ).map{ |ln| ln.split('; ')[1] }
      end
      genomes.uniq
   end
   def genome2taxid(genome_id)
      ln = ebiFetch('embl', [genome_id], 'annot').split(/[\n\r]/).grep(/^FT\s+\/db_xref="taxon:/).first
      return ln if ln.nil?
      ln.sub(/.*"taxon:(\d+)".*/, "\\1")
   end
   def genome2taxon(genome_id, rank='species')
      xml = ebiFetch('taxonomy', [genome2taxid(genome_id)], 'enataxonomyxml').gsub(/\s*\n\s*/,'')
      xml.scan(/<taxon [^>]+>/).grep(/rank="#{rank}"/).first.sub(/.* taxId="(\d+)".*/,"\\1")
   end
   def restcall(url, outfile=nil)
      response = RestClient::Request.execute(:method=>:get,  :url=>url, :timeout=>600)
      raise "Unable to reach EBI REST client, error code #{response.code}." unless response.code == 200
      unless outfile.nil?
	 ohf = File.open(outfile, 'w')
	 ohf.print response.to_s
	 ohf.close
      end
      response.to_s
   end
   def ebiFetch(db, ids, format, outfile=nil)
      url = "#{ROCker.ebirest}/dbfetch/dbfetch/#{db.to_s}/#{ids.join(",")}/#{format.to_s}"
      $stderr.puts "   # Calling: #{url}" if @o[:debug]
      res = self.restcall url
      unless outfile.nil?
	 ohf = File.open(outfile, 'w')
	 ohf.print res
	 ohf.close
      end
      res
   end
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

