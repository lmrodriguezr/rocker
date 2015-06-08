#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jun-08-2015
#

class ROCker
   #================================[ Class ]
   #@@DEFAULTS.merge!({  })
   
   #================================[ Filter ]
   def filter!(data=nil)
      raise "-k/--rocker is mandatory." if @o[:rocker].nil?
      raise "-x/--query-blast is mandatory." if @o[:qblast].nil?
      raise "-o/--out-blast is mandatory." if @o[:oblast].nil?
      
      # Read ROCker file
      if data.nil?
	 puts "Loading ROCker file: #{@o[:rocker]}." unless @o[:q]
	 data = ROCData.new @o[:rocker]
      end

      # Filter similarity search
      puts "Filtering similarity search: #{@o[:qblast]}." unless @o[:q]
      ih = File.open(@o[:qblast], "r")
      oh = File.open(@o[:oblast], "w")
      while ln = ih.gets
	 bh = BlastHit.new(ln, data.aln)
	 oh.print ln if not(bh.sfrom.nil?) and bh.bits >= data.win_at_col(bh.midpoint).thr
      end
      ih.close
      oh.close
   end # filter!
end # ROCker

