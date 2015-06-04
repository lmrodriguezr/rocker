#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jun-04-2015
#

class ROCker
   #================================[ Class ]
   #@@DEFAULTS.merge!({  })
   
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
end # ROCker

