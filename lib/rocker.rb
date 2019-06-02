#
# @author  Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author  Luis (Coto) Orellana
# @license Artistic-2.0
#

require "rocker/blasthit"
require "rocker/rocdata"

class ROCker
   #================================[ Class ]
   @@VERSION = "1.3.1"
   @@CITATION = [
      "Orellana, Rodriguez-R & Konstantinidis, 2016. DOI:10.1093/nar/gkw900.",
      "ROCker: accurate detection and quantification of target genes in",
      "short-read metagenomic data sets by modeling sliding-window bitscores.",
      "Nucleic Acids Research 45(3):e14."]
   @@DATE = "2019-06-02"
   @@DEFAULTS = {
      # General
      q: false, r: "R", nucl: false, debug: false, thr: 2, search: :blast,
      # External software
      searchbins: "",
      searchcmd: {
	 blast: '%1$s%2$s -query "%3$s" -db "%4$s" -out "%5$s" ' +
	    '-num_threads %6$d -outfmt 6 -max_target_seqs 1',
	 diamond: '%1$sdiamond %2$s -q "%3$s" -d "%4$s" -a "%5$s.daa" -p %6$d' +
	    ' -k 1 --min-score 20 --sensitive && %1$sdiamond view -a "%5$s"' +
	    ' -o "%5$s"'},
      makedbcmd: {
	 blast: '%1$smakeblastdb -dbtype %2$s -in "%3$s" -out "%4$s"',
	 diamond: '%1$sdiamond makedb --in "%3$s" -d "%4$s"'}
   }
   def self.defaults() @@DEFAULTS ; end
   def self.default(k) @@DEFAULTS[k] ; end
   def self.VERSION; @@VERSION ; end
   def self.DATE; @@DATE ; end
   def self.CITATION(j=" ") @@CITATION.join(j) ; end

   #================================[ Instance ]
   attr_reader :o
   def initialize(opts)
      @o = ROCker.defaults
      opts.each{ |k,v| @o[k] = v }
      RInterface.R_BIN = opts[:r] unless opts[:r].nil?
   end
   
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
   def bash(cmd, err_msg=nil)
      o = `#{cmd} 2>&1 && echo '{'`
      raise (err_msg.nil? ? "Error executing: #{cmd}\n\n#{o}" : err_msg) unless
	 o[-2]=="{"
      true
   end
end

#================================[ Extensions ]
# To ROCker
require "rocker/step/build"
require "rocker/step/compile"
require "rocker/step/search"
require "rocker/step/filter"
require "rocker/step/plot"

# To other
class Numeric
   def ordinalize
      n= self.to_s
      s= n[-2]=='1' ? "th" :
	 n[-1]=='1' ? "st" :
	 n[-1]=='2' ? "nd" :
	 n[-1]=='3' ? "rd" : "th"
      n + s
   end
end

