#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Dec-01-2015
#

require "rocker/alignment"

class ProteinSet
   attr_reader :rocker, :ids, :aln
   def initialize(rocker, ids=nil, file=nil, aln_file=nil)
      @genomes = {}
      @tranids = {}
      @aln = nil
      @rocker = rocker
      @ids = []
      @ids += ids unless ids.nil?
      @ids += File.readlines(file).map{ |l| l.chomp } unless file.nil?
      unless aln_file.nil?
	 aln = Alignment.new
	 aln.read_fasta aln_file
	 aln_ids = aln.get_ids
	 @aln = aln if (@ids - aln_ids).empty?
	 @ids += aln_ids
      end
      @ids.uniq!
   end
   def download(file)
      tmp_ids = Array.new(self.ids)
      f = File.open(file, "w")
      while tmp_ids.size>0
	 f.print rocker.ebiFetch(:uniprotkb, tmp_ids.shift(200), :fasta)
      end
      f.close
   end
   def get_from_aln(file, aln)
      f = File.open(file, "w")
      f.print aln.to_seq_s
      f.close
   end
   def get_genomes!
      self.ids.each do |id|
	 doc = self.rocker.ebiFetch(:uniprotkb, [id], :annot).split("\n")
	 doc.grep( /^DR\s+EMBL;/ ).map do |ln|
	    r=ln.split("; ")
	    self.link_genome(id, r[1])
	    self.link_tranid(id, r[2])
	 end
      end
   end
   def link_genome(prot_id, genome_id)
      @genomes[prot_id] ||= []
      @genomes[prot_id] << genome_id
      @genomes[prot_id].uniq!
   end
   def link_tranid(prot_id, transl_id)
      @tranids[prot_id] ||= []
      @tranids[prot_id] << transl_id
      @tranids[prot_id].uniq!
   end
   def genomes
      return [] if @genomes.empty?
      @genomes.values.reduce(:+).uniq
   end
   def tranids
      return [] if @tranids.empty?
      @tranids.values.reduce(:+).uniq
   end
   def tranids_dump
      @tranids.map{|k,v| "{#{k}: #{v}}"}.join(", ")
   end
   def in_coords(coords)
      coords.keys.map do |genome|
	 locations = coords[ genome ]
	 locations.map do |loc|
	    if not loc[:prot_id].nil? 
	       loc[:prot_id] if include? loc[:prot_id]
	    elsif not loc[:tran_id].nil?
	       @tranids.map{ |k,v| v.include?(loc[:tran_id]) ? k : nil }.compact.first
	    else
	       warn "Warning: Impossible to resolve protein located " +
		  "in '#{genome}' at: #{loc}."
	       nil
	    end
	 end
      end.reduce([], :+).compact.uniq
   end
   def size() self.ids.size end
   def empty?() self.ids.empty? end
   def include?(id) self.ids.include?(id) end
end


