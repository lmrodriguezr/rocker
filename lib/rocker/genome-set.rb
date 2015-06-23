#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jun-23-2015
#

class GenomeSet
   attr_reader :rocker, :ids, :taxa
   def initialize(rocker, ids)
      @rocker = rocker
      @ids = ids
      @ids = [] if ids.nil?
      @taxa = {}
      @all_taxa = {}
   end
   def download(file)
      tmp_ids = Array.new(self.ids)
      ofh = File.open(file, "w")
      while tmp_ids.size>0
	 ofh.print rocker.ebiFetch(:embl, tmp_ids.shift(200), :fasta)
      end
      ofh.close
   end
   def link_taxon(id, taxon)
      @all_taxa[ taxon.to_sym ] ||= []
      @all_taxa[ taxon.to_sym ] << id
   end
   def choose_genomes!(rank)
      @taxa = {}
      self.get_taxonomy! rank
      @all_taxa.each_pair{ |taxon,ids| @taxa[taxon] = ids.sample }
      @ids = @taxa.values
   end
   def get_taxonomy!(rank)
      @all_taxa = {}
      ids.each do |id|
	 self.link_taxon(id, genome2taxon(id, rank))
      end
   end
   def taxa=(hash)
      @taxa = {}
      hash.each_pair{ |taxon, id| @taxa[taxon] = id if self.ids.include? id }
   end
   def size() self.ids.size end
   def empty?() self.ids.empty? end

   #================================[ Utilities ]
   def genome2taxon(genome_id, rank='species')
      v = genome2taxid(genome_id)
      unless v.nil?
	 xml = rocker.ebiFetch('taxonomy', [v], 'enataxonomyxml').gsub(/\s*\n\s*/,'')
	 v = xml.scan(/<taxon [^>]+>/).grep(/rank="#{rank}"/).first
	 v.sub!(/.* taxId="(\d+)".*/,"\\1") unless v.nil?
      end
      return "no-taxon-#{(0...12).map { (65 + rand(26)).chr }.join}" if v.nil? or v !~ /^\d+$/
      v
   end
   def genome2taxid(genome_id)
      doc = rocker.ebiFetch('embl', [genome_id], 'annot').split(/[\n\r]/)
      ln = doc.grep(/^FT\s+\/db_xref="taxon:/).first
      ln = doc.grep(/^OX\s+NCBI_TaxID=/).first if ln.nil?
      return nil if ln.nil?
      ln.sub!(/.*(?:"taxon:|NCBI_TaxID=)(\d+)["; ].*/, "\\1")
      return nil unless ln =~ /^\d+$/
      ln
   end
end


