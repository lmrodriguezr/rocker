#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jun-04-2015
#

require 'json'

class ROCker
   #================================[ Class ]
   @@EBIREST = 'http://www.ebi.ac.uk/Tools'
   @@DEFAULTS.merge!({:positive=>[], :negative=>[], :genomefrx=>1.0, :seqdepth=>0.03, :readlen=>100, :minovl=>50,
      # Ext. Software
      :grinder=>'grinder', :muscle=>'muscle',
      :grindercmd=>'%1$s -reference_file "%2$s" -cf "%3$f" -dc \'-~*NnKkMmRrYySsWwBbVvHhDdXx\' -md uniform 0.1 -mr 95 5 -rd %4$d uniform 5 -base_name "%5$s"',
      :musclecmd=>'%1$s -in "%2$s" -out "%3$s" -quiet'
   })
   @@HAS_BUILD_GEMS = nil
   def self.ebirest() @@EBIREST ; end
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

   #================================[ Utilities ]
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
      $stderr.puts "   # Calling: #{url}" if @o[:debug]
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
      res = self.restcall url
      unless outfile.nil?
	 ohf = File.open(outfile, 'w')
	 ohf.print res
	 ohf.close
      end
      res
   end
   def get_coords_from_gff3(genome_ids, protein_ids, thread_id)
      positive_coords = {}
      genomes_org = {}
      i = 0
      genome_ids.each do |genome_id|
	 print "  * scanning #{(i+=1).ordinalize} genome out of #{genome_ids.size} [thread #{thread_id}]. Threads done: #{@o[:threads_done]}.     \r" unless @o[:q]
	 unless @o[:pertaxon].nil?
	    genome_taxon = genome2taxon(genome_id, @o[:pertaxon])
	    genomes_org[ genome_taxon.to_sym ] ||= []
	    genomes_org[ genome_taxon.to_sym ] << genome_id
	 end
	 genome_file = @o[:baseout] + ".src." + genome_id + ".gff3"
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
	    positive_coords[ r[0].to_sym ] ||= []
	    positive_coords[ r[0].to_sym ] << {
	       :prot_id	=> p,
	       :from	=> r[3].to_i,
	       :to	=> r[4].to_i,
	       :strand	=> r[6]
	    }
	 end
      end
      {:positive_coords=>positive_coords, :genomes_org=>genomes_org}
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
	    puts "  * linking genomes from #{@o[set].size} #{set.to_s} sequence(s)." unless @o[:q]
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
	 thrs = [@o[:thr], genome_ids[:positive].size].min
	 puts "  * downloading and parsing #{genome_ids[:positive].size} GFF3 document(s) in #{thrs} threads." unless @o[:q]
	 $stderr.puts "   # Looking for: #{@o[:positive]}" if @o[:debug]
	 $stderr.puts "   # Looking into: #{genome_ids[:positive]}" if @o[:debug]
	 thr_obj = []
	 @o[:threads_done] = 0
	 (0 .. (thrs-1)).each do |thr_i|
	    thr_obj << Thread.new do
	       Thread.current[:ids_to_parse] = []
	       Thread.current[:i] = 0
	       while Thread.current[:i] < genome_ids[:positive].size
		  Thread.current[:ids_to_parse] << genome_ids[:positive][ Thread.current[:i] ] if (Thread.current[:i] % thrs)==thr_i
		  Thread.current[:i] += 1
	       end
	       Thread.current[:output] = get_coords_from_gff3(Thread.current[:ids_to_parse], genome_ids[:positive], thr_i)
	       @o[:threads_done] += 1
	    end
	 end
	 # Combine results
	 positive_coords = {}
	 genomes_org = {}
	 genome_org = {}
	 thr_obj.each do |t|
	    t.join
	    raise "Thread failed without error trace: #{t}" if t[:output].nil?
	    t[:output][:positive_coords].each_pair do |k,v|
	       positive_coords[ k ] ||= []
	       positive_coords[ k ] += v
	    end
	    t[:output][:genomes_org].each_pair do |k,v|
	       genomes_org[ k ] ||= []
	       genomes_org[ k ] << v
	    end
	 end
	 print "\n" unless @o[:q]

	 # Select one genome per taxon
	 unless @o[:pertaxon].nil?
	    genomes_org.each_pair{ |k,v| genome_org[ k ] = v.sample.first }
	 end
	 
	 # Save coordinates
	 ofh = File.open(coords_file, "w")
	 ofh.print JSON.pretty_generate({:positive_coords=>positive_coords, :genome_org=>genome_org})
	 ofh.close
      end
      unless @o[:pertaxon].nil?
	 genome_ids[:positive] = genome_org.values
	 puts "  Using #{genome_org.size} genome(s) after filtering by #{@o[:pertaxon]}." unless @o[:q]
      end
      all_genome_ids = genome_ids.values.reduce(:+).uniq
      found = positive_coords.values.map{ |a| a.map{ |b| b[:prot_id] } }.reduce(:+)
      raise "Cannot find the genomic location of any provided sequence." if found.nil?
      missing = @o[:positive] - found
      warn "\nWARNING: Cannot find genomic location of sequence(s) #{missing.join(',')}.\n\n" unless missing.size==0 or @o[:genomefrx]<1.0
      
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
	    puts "  * running grinder and tagging positive reads in #{thrs} threads." unless @o[:q]
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
			raise "Cannot parse simulated read's defline, are you using Grinder?: #{Thread.current[:l]}" if Thread.current[:rd].nil?
			Thread.current[:positive] = false
			positive_coords[Thread.current[:rd][:genome_id].to_sym] ||= []
			positive_coords[Thread.current[:rd][:genome_id].to_sym].each do |gn|
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
end # ROCker

