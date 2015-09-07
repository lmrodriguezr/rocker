#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Sep-07-2015
#

class ROCWindow
   attr_reader :data, :from, :to, :hits, :tps, :thr
   def initialize(data, from=nil, to=nil)
      @data = data
      if from.is_a? String
	 r = from.split(/\t/)
	 @from	= r[0].to_i
	 @to	= r[1].to_i
	 @hits	= r[2].to_i
	 @tps	= r[3].to_i
	 @thr	= r[4].to_f
      else
	 a = from.nil? ? 1 : [from,1].max
	 b = to.nil? ? data.aln.cols : [to,data.aln.cols].min
	 @from = [a,b].min
	 @to = [a,b].max
	 @thr = nil
	 compute!
      end
   end
   def compute!
      load_hits
      @hits = rrun("nrow(y);", :int)
      @tps = rrun("sum(y$V5==1);", :int)
      unless almost_empty
	 rrun "rocobj <- roc(as.numeric(y$V5==1), y$V4);"
	 thr = rrun("coords(rocobj, 'best', ret='threshold', " +
	    "best.method='youden', " +
	    "best.weights=c(0.5, sum(y$V5==1)/nrow(y)))[1];", :float)
	 @thr = thr.to_f
	 @thr = nil if @thr==0.0 or @thr.infinite?
      end
   end
   def around_thr
      a = self.previous
      b = self.next
      while not a.nil? and a.thr.nil?
	 a = a.previous
      end
      while not b.nil? and b.thr.nil?
	 b = b.next
      end
      return nil if a.nil? and b.nil?
      return a.thr if b.nil?
      return b.thr if a.nil?
      return (b.thr*(from-a.from) - a.thr*(from-b.from))/(b.from-a.from)
   end
   def load_hits() self.rrun "y <- x[x$V6>=#{from} & x$V6<=#{to},];" end
   def previous() (from == 1) ? nil : data.win_at_col(from - 1) end
   def next() (to == data.aln.cols) ? nil : data.win_at_col(to + 1) end
   def thr_notnil() (@thr.nil? or @thr.infinite?) ? around_thr : @thr end
   def fps() hits - tps end
   def almost_empty() fps < 3 or tps < 3 end
   def length() to - from + 1 end
   def rrun(cmd, type=nil) data.rrun(cmd, type) end
   def to_s() [from, to, hits, tps, thr_notnil].join("\t") + "\n" end
end

