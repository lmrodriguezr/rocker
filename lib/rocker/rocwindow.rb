#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jan-22-2015
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
	 self.compute!
      end
   end
   def compute!
      self.load_hits
      @hits = self.rrun "nrow(y);", :int
      @tps = self.rrun "sum(y$V5);", :int
      unless self.almost_empty
	 self.rrun "rocobj <- roc(y$V5, y$V4);"
	 thr = self.rrun 'coords(rocobj, "best", ret="threshold", best.method="youden", best.weights=c(0.5, sum(y$V5)/nrow(y)))[1];', :float
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
      return (b.thr*(self.from-a.from) - a.thr*(self.from-b.from))/(b.from-a.from)
   end
   def load_hits() self.rrun "y <- x[x$V6>=#{self.from} & x$V6<=#{self.to},];" end
   def previous() (self.from == 1) ? nil : self.data.win_at_col(self.from - 1) end
   def next() (self.to == self.data.aln.cols) ? nil : self.data.win_at_col(self.to + 1) end
   def thr_notnil() (@thr.nil? or @thr.infinite?) ? self.around_thr : @thr end
   def fps() self.hits - self.tps end
   def almost_empty() self.fps < 3 or self.tps < 3 end
   def length() self.to - self.from + 1 end
   def rrun(cmd, type=nil) self.data.rrun cmd, type end
   def to_s() [self.from, self.to, self.hits, self.tps, self.thr_notnil].join("\t") + "\n" end
end

