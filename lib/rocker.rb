#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jun-05-2015
#

require 'rocker/blasthit'
require 'rocker/rocdata'

class ROCker
   #================================[ Class ]
   @@DEFAULTS = {
      # General
      :q=>false, :r=>'R', :nucl=>false, :debug=>false,:thr=>2,:search=>:blast,
      # External software
      :searchbins=>'',
      :searchcmd=>{
	 :blast=>'%1$s%2$s -query "%3$s" -db "%4$s" -out "%5$s" -num_threads %6$d -outfmt 6 -max_target_seqs 1',
	 :diamond=>'%1$sdiamond %2$s -q "%3$s" -d "%4$s" -o "%5$s" -t %6$d -k 1 --min-score 20 --sensitive'},
      :makedbcmd=>{
	 :blast=>'%1$smakeblastdb -dbtype %2$s -in "%3$s" -out "%4$s"',
	 :diamond=>'%1$sdiamond makedb --in "%3$s" -d "%4$s"'}
   }
   def self.defaults() @@DEFAULTS ; end
   def self.default(k) @@DEFAULTS[k] ; end

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
      raise (err_msg.nil? ? "Error executing: #{cmd}\n\n#{o}" : err_msg) unless o[-2]=='{'
      true
   end
end

#================================[ Extensions ]
# To ROCker
require 'rocker/step/build'
require 'rocker/step/compile'
require 'rocker/step/search'
require 'rocker/step/filter'
require 'rocker/step/plot'

# To other
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

