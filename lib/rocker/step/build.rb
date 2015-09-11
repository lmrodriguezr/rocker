#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Sep-11-2015
#

require 'json'
require 'rocker/protein-set'
require 'rocker/genome-set'

class ROCker
   #================================[ Class ]
   @@EBIREST = "http://www.ebi.ac.uk/Tools"
   @@DEFAULTS.merge!({positive:[], negative:[], seqdepth:0.03, readlen:100,
      minovl:50,
      # Ext. Software
      aligner: :clustalo, simulator: :grinder,
      simulatorbin:{grinder:"grinder"},
      simulatorcmd:{grinder:"%1$s -reference_file \"%2$s\" -cf \"%3$f\" " +
	 "-dc '-~*NnKkMmRrYySsWwBbVvHhDdXx' -md uniform 0.1 -mr 95 5 " +
	 "-rd %4$d uniform 5 -base_name \"%5$s\""},
      alignerbin:{muscle:"muscle", clustalo:"clustalo"},
      alignercmd:{muscle:"%1$s -in \"%2$s\" -out \"%3$s\" -quiet",
	 clustalo:"%1$s -i \"%2$s\" -o \"%3$s\" --threads=%4$d --force"}
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
   def restcall(url, outfile=nil)
      $stderr.puts "   # Calling: #{url}" if @o[:debug]
      response = RestClient::Request.execute(:method=>:get, :url=>url,
	 :timeout=>600)
      raise "Unable to reach EBI REST client, error code " +
	 response.code.to_s + "." unless response.code == 200
      unless outfile.nil?
	 ohf = File.open(outfile, "w")
	 ohf.print response.to_s
	 ohf.close
      end
      response.to_s
   end
   def ebiFetch(db, ids, format, outfile=nil)
      url = "#{ROCker.ebirest}/dbfetch/dbfetch/" +
	 "#{db.to_s}/#{ids.join(",")}/#{format.to_s}"
      self.restcall url, outfile
   end
   def get_coords_from_gff3(genome_ids, pset, thread_id, json_file)
      coords = {}
      i = 0
      genome_ids.each do |genome_id|
	 print "  * scanning #{(i+=1).ordinalize} genome out of " +
	    "#{genome_ids.size} in first thread.  \r" if
	    thread_id==0 and not @o[:q]
	 genome_file = @o[:baseout] + ".src." + genome_id + ".gff3"
	 if @o[:reuse] and File.size? genome_file
	    ifh = File.open(genome_file, "r")
	    doc = ifh.readlines.grep(/^[^#]/)
	    ifh.close
	 else
	    genome_file=nil unless @o[:noclean]
	    doc = ebiFetch(:embl, [genome_id], :gff3,
	       genome_file).split("\n").grep(/^[^#]/)
	 end
	 doc.each do |ln|
	    next if ln =~ /^#/
	    r = ln.chomp.split /\t/
	    next if r.size < 9
	    prots = r[8].split(/;/).grep(
	       /^db_xref=UniProtKB[\/A-Za-z-]*:/){ |xref| xref.split(/:/)[1] }
	    p = prots.select{ |id| pset.ids.include? id }.first
	    trans = r[8].split(/;/).grep(
	       /^protein_id=/){ |pid| pid.split(/=/)[1] }
	    t = trans.select{ |id| pset.tranids.include? id }.first
	    next if p.nil? and t.nil?
	    coords[ r[0].to_sym ] ||= []
	    coords[ r[0].to_sym ] << {
	       prot_id:	p,
	       tran_id:	t,
	       from:	r[3].to_i,
	       to:	r[4].to_i,
	       strand:	r[6]
	    }
	 end
      end
      print "\n" if thread_id==0 and not @o[:q]
      ofh = File.open(json_file, "w")
      ofh.print({coords:coords}.to_json)
      ofh.close
   end
   
   #================================[ Build ]
   def build!
      # Check requirements
      puts "Testing environment." unless @o[:q]
      {  searchcmd: :search, makedbcmd: :search,
	 alignercmd: :aligner, alignerbin: :aligner,
	 simulatorcmd: :simulator, simulatorbin: :simulator
      }.each_pair { |k,v| @o[k] = @o[k][@o[v]] if @o[k].is_a? Hash }
      @o[:nosearch]=true if @o[:nosimulate]
      raise "Unsatisfied requirements, please see the help message (-h)." unless
	 ROCker.has_build_gems?
      protein_set = {}
      protein_set[:+] = ProteinSet.new(self,@o[:positive],@o[:posfile],@o[:aln])
      protein_set[:-] = ProteinSet.new(self,@o[:negative],@o[:negfile])
      raise "-p, -P, or -a are mandatory." if protein_set[:+].empty?
      raise "-o/--baseout is mandatory." if @o[:baseout].nil?
      if protein_set[:+].size==1 and not @o[:noaln]
	 warn "\nWARNING: Positive set contains only one sequence, turning " +
	    "off alignment.\n\n"
	 @o[:noaln] = true
      end
      unless @o[:nosimulate]
	 self.bash("#{@o[:simulatorbin]} --version",
	    "--simulator-bin must be executable. Is Grinder installed?") if
	    @o[:simulator]==:grinder
      end
      unless @o[:noaln]
	 self.bash("#{@o[:alignerbin]} -version",
	    "--aligner-bin must be executable. Is Muscle installed?") if
	    @o[:aligner]==:muscle
	 self.bash("#{@o[:alignerbin]} --version",
	    "--aligner-bin must be executable. Is ClustalOmega installed?") if
	    @o[:aligner]==:clustalo
      end
      unless @o[:nosearch]
	 self.bash("#{@o[:searchbins]}makeblastdb -version",
	    "--search-bins must contain executables. Is BLAST+ installed?") if
	    @o[:search]==:blast
	 self.bash("#{@o[:searchbins]}diamond --help",
	    "--search-bins must contain executables. Is DIAMOND installed?") if
	    @o[:search]==:diamond
      end

      # Download genes
      puts "Downloading gene data." unless @o[:q]
      ref_file = @o[:baseout] + ".ref.fasta"
      if not protein_set[:+].aln.nil?
	 puts "  * reusing aligned sequences as positive set." unless @o[:q]
	 protein_set[:+].get_from_aln(ref_file, aln)
	 @o[:noaln] = true
      elsif @o[:reuse] and File.size? ref_file
	 puts "  * reusing positive set: #{ref_file}." unless @o[:q]
      else
	 puts "  * downloading #{protein_set[:+].size} sequence(s) in " +
	    "positive set." unless @o[:q]
	 $stderr.puts "   # #{protein_set[:+].ids}" if @o[:debug]
	 protein_set[:+].download(ref_file)
      end
      [:+, :-].each do |set|
         unless protein_set[set].empty?
	    puts "  * linking genomes from #{protein_set[set].size} " +
	       "[#{set.to_s}] sequence(s)." unless @o[:q]
	    $stderr.puts "   # #{protein_set[set].ids}" if @o[:debug]
	    protein_set[set].get_genomes!
	 end
      end
      raise "No genomes associated with the positive set." if
	 protein_set[:+].genomes.empty?
      genome_set = {:+ => GenomeSet.new(self, protein_set[:+].genomes),
	 :- => GenomeSet.new(self, protein_set[:-].genomes)}
      
      # Locate genes
      puts "Analyzing genome data." unless @o[:q]
      coords_file = @o[:baseout] + ".src.coords"
      if @o[:reuse] and File.size? coords_file
	 puts "  * reusing coordinates: #{coords_file}." unless @o[:q]
	 c = JSON.parse File.read(coords_file), {symbolize_names:true}
	 positive_coords = c[:positive_coords]
	 negative_coords = c[:negative_coords]
	 genome_set[:+].taxa = c[:taxa_pos]
	 genome_set[:-].taxa = c[:taxa_neg]
      else
	 all_coords = {}
	 [:+, :-].each do |set_type|
	    all_coords[set_type] = {}
	    next if genome_set[set_type].empty?
	    thrs = [@o[:thr], genome_set[set_type].size].min
	    puts "  * downloading and parsing #{genome_set[set_type].size} " +
	       "GFF3 document(s) in #{thrs} threads." unless @o[:q]
	    $stderr.puts "   # Looking for translations: " +
	       "#{protein_set[set_type].tranids}" if @o[:debug]
	    $stderr.puts "   # Looking into: #{genome_set[set_type].ids}" if
	       @o[:debug]
	    # Launch threads
	    thr_obj = []
	    (0 .. (thrs-1)).each do |thr_i|
	       ids_to_parse = []
	       (0 .. (genome_set[set_type].size-1)).each do |i|
		  ids_to_parse << protein_set[set_type].genomes[i] if
		     (i % thrs) == thr_i
	       end
	       json_file = @o[:baseout] + ".src.coords." + thr_i.to_s + ".tmp"
	       thr_obj << json_file
	       fork do
		  get_coords_from_gff3(ids_to_parse, protein_set[set_type],
		     thr_i, json_file)
	       end
	    end
	    # Combine results
	    Process.waitall
	    thr_obj.each do |t|
	       raise "Thread failed without error trace: #{t}" unless
		  File.exist? t
	       o = JSON.parse(File.read(t), {symbolize_names:true})
	       o[:coords].each_pair do |k,v|
		  all_coords[set_type][ k ] ||= []
		  all_coords[set_type][ k ] += v
	       end
	       File.unlink t
	    end
	 end # [:+, :-].each
	 positive_coords = all_coords[:+]
	 negative_coords = all_coords[:-]
	 # Select one genome per taxon
	 unless @o[:pertaxon].nil?
	    puts "  Selecting genomes by #{@o[:pertaxon]}." unless @o[:q]
	    [:+,:-].each{ |set| genome_set[set].choose_genomes! @o[:pertaxon] }
	 end
	 # Save coordinates and taxa
	 ofh = File.open(coords_file, "w")
	 ofh.print JSON.pretty_generate({
	    positive_coords:positive_coords,
	    negative_coords:negative_coords,
	    taxa_pos:genome_set[:+].taxa,
	    taxa_neg:genome_set[:-].taxa})
	 ofh.close
      end # if @o[:reuse] and File.size? coords_file ... else
      unless @o[:pertaxon].nil?
	 puts "  Using " +
	    [:+,:-].map{ |set| genome_set[set].size }.reduce(:+).to_s +
	    " genome(s) after filtering by #{@o[:pertaxon]}." unless @o[:q]
      end
      found = protein_set[:+].in_coords(positive_coords)
      raise "Cannot find the genomic location of any provided sequence." if
	 found.nil?
      missing = protein_set[:+].ids - found
      warn "\nWARNING: Cannot find genomic location of #{missing.size} " +
	 "sequence(s) #{missing.join(",")}.\n\n" unless missing.empty?
      
      # Download genomes
      genome_set[:all] = GenomeSet.new(self,
	 genome_set[ :+ ].ids + genome_set[ :- ].ids)
      genomes_file = @o[:baseout] + ".src.fasta"
      if @o[:reuse] and File.size? genomes_file
	 puts "  * reusing existing file: #{genomes_file}." unless @o[:q]
      else
	 puts "  * downloading " + genome_set[:all].size.to_s +
	    " genome(s) in FastA." unless @o[:q]
	 $stderr.puts "   # #{genome_set[:all].ids}" if @o[:debug]
	 genome_set[:all].download genomes_file
      end
      
      # Generate metagenome
      unless @o[:nosimulate]
	 puts "Generating in silico metagenome" unless @o[:q]
	 if @o[:reuse] and File.size? @o[:baseout] + ".mg.fasta"
	    puts "  * reusing existing file: #{@o[:baseout]}.mg.fasta." unless
	       @o[:q]
	 else
	    all_src = File.readlines("#{@o[:baseout]}.src.fasta"
	       ).select{ |l| l =~ /^>/ }.size
	    thrs = [@o[:thr], all_src].min
	    thr_obj = []
	    seqs_per_thr = (all_src.to_f/thrs).ceil
	    thrs = (all_src.to_f/seqs_per_thr).ceil
	    puts "  * simulating metagenomes and tagging positive reads in " +
	       thrs.to_s + " threads." unless @o[:q]
	    $stderr.puts "   # #{positive_coords}" if @o[:debug]
	    (0 .. (thrs-1)).each do |thr_i|
	       output = @o[:baseout] + ".mg.fasta.#{thr_i.to_s}"
	       thr_obj << output
	       fork do
		  seqs_a = thr_i*seqs_per_thr + 1
		  seqs_b = [seqs_a + seqs_per_thr - 1, all_src].min
		  # Create sub-fasta
		  ofh = File.open("#{@o[:baseout]}.src.fasta.#{thr_i.to_s}","w")
		  ifh = File.open("#{@o[:baseout]}.src.fasta","r")
		  seq_i = 0
		  while l = ifh.gets
		     seq_i+=1 if l =~ /^>/
			     break if seq_i > seqs_b
		     ofh.print l if seq_i >= seqs_a
		  end
		  ifh.close
		  ofh.close

		  # Run simulator (except if the temporal file is already
		  # there and can be reused)
		  bash sprintf(@o[:simulatorcmd], @o[:simulatorbin],
		     "#{@o[:baseout]}.src.fasta.#{thr_i.to_s}",
		     @o[:seqdepth]*@o[:readlen].to_f, @o[:readlen],
		     "#{@o[:baseout]}.mg.tmp.#{thr_i.to_s}") unless
			@o[:reuse] and
			File.size? @o[:baseout] +
			".mg.tmp.#{thr_i.to_s}-reads.fa"

		  # Tag positive and negative reads
		  puts "  * tagging reads [thread #{thr_i}]." unless
		     @o[:q]
		  ifh = File.open(@o[:baseout] + ".mg.tmp.#{thr_i}-reads.fa",
		     "r")
		  ofh = File.open(@o[:baseout] + ".mg.fasta.#{thr_i}", "w")
		  while l = ifh.gets
		     if l =~ /^>/
			rd = %r{
			   ^>(?<id>\d+)\s
			   reference=[A-Za-z]+\|
			   (?<genome_id>[A-Za-z0-9_]+)\|.*\s
			   position=(?<comp>complement\()?(?<from>\d+)\.\.
			   (?<to>\d+)\)?\s
			}x.match(l)
			raise "Cannot parse simulated read's defline, are " +
			   "you using Grinder?: #{l}" if rd.nil?
			positive = false
			positive_coords[rd[:genome_id].to_sym] ||= []
			positive_coords[rd[:genome_id].to_sym].each do |gn|
			   left  = rd[:to].to_i - gn[:from]
			   right = gn[:to] - rd[:from].to_i
			   if (left*right >= 0) and
				 ([left, right].min >= @o[:minovl])
			      positive = true
			      break
			   end
			end
			negative = false
			negative_coords[rd[:genome_id].to_sym] ||= []
			negative_coords[rd[:genome_id].to_sym].each do |gn|
			   left  = rd[:to].to_i - gn[:from]
			   right = gn[:to] - rd[:from].to_i
			   if (left*right >= 0) and
				 ([left, right].min >= @o[:minovl])
			      negative = true
			      break
			   end
			end
			l = ">#{thr_i.to_s}_#{rd[:id]}" +
			   "#{positive ? "@%" : (negative ? "@$" : "")} " +
			   "ref=#{rd[:genome_id]}:#{rd[:from]}..#{rd[:to]}" +
			   "#{(rd[:comp]=="complement(") ? "-" : "+"}\n"
		     end
		     ofh.print l
		  end
		  ofh.close
		  ifh.close
	       end # fork
	    end # (1 .. thrs).each
	    Process.waitall
	    # Concatenate results
	    ofh = File.open(@o[:baseout] + ".mg.fasta", "w")
	    thr_obj.each do |t|
	       raise "Thread failed without error trace: #{t}" unless
		  File.exist? t
	       ifh = File.open(t, "r")
	       while l = ifh.gets
	          ofh.print l
	       end
	       ifh.close
	       File.unlink t
	    end
	    ofh.close
         end
      end # unless @o[:nosimulate]
      
      # Align references
      unless @o[:noaln]
	 puts "Aligning reference set." unless @o[:q]
	 if @o[:reuse] and File.size? "#{@o[:baseout]}.ref.aln"
	    puts "  * reusing existing file: #{@o[:baseout]}.ref.aln." unless
	       @o[:q]
	 else
	    bash(sprintf(@o[:alignercmd],
	       @o[:alignerbin], "#{@o[:baseout]}.ref.fasta",
	       "#{@o[:baseout]}.ref.aln", @o[:thr]))
	    puts "  +--\n  | IMPORTANT NOTE: Manually checking the alignment " +
	       "before\n  | the 'compile' step is *strongly* encouraged.\n  " +
	       "+--\n" unless @o[:q]
	 end
      end
      
      # Run similarity search
      unless @o[:nosearch]
	 puts "Running similarity search." unless @o[:q]
	 if @o[:reuse] and File.size? "#{@o[:baseout]}.ref.blast"
	    puts "  * reusing existing file: #{@o[:baseout]}.ref.blast." unless
	       @o[:q]
	 else
	    puts "  * preparing database." unless @o[:q]
	    bash(sprintf(@o[:makedbcmd],
	       @o[:searchbins], "prot", "#{@o[:baseout]}.ref.fasta",
	       "#{@o[:baseout]}.ref"))
	    puts "  * running similarity search." unless @o[:q]
	    bash(sprintf(@o[:searchcmd],
	       @o[:searchbins], "blastx", "#{@o[:baseout]}.mg.fasta",
	       "#{@o[:baseout]}.ref", "#{@o[:baseout]}.ref.blast", @o[:thr]))
	 end
      end
      
      # Clean
      unless @o[:noclean]
	 puts "Cleaning." unless @o[:q]
	 sff  = %w{.src.xml .src.fasta}
	 sff += %w{.mg.tmp-reads.fa .mg.tmp-ranks.txt} unless @o[:nosimulate]
	 sff += %w{.ref.phr .ref.pin .ref.psq} unless @o[:nosearch]
	 sff.each { |sf| File.unlink @o[:baseout] + sf if
	    File.exist? @o[:baseout] + sf }
      end
   end # build!
end # ROCker

