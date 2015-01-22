#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jan-22-2015
#

require 'rocker/rinterface'
require 'rocker/rocwindow'
require 'rocker/alignment'
require 'tmpdir'

class ROCData
   attr_reader :aln, :windows, :r
   # Use ROCData.new(table,aln,window) to re-compute from table, use ROCData.new(data) to load
   def initialize(val, aln=nil, window=nil)
      @r = RInterface.new
      @nucl = false
      if not aln.nil?
	 @aln = aln
	 self.rrun "library('pROC');"
	 self.rrun "x <- read.table('#{val}', sep='\\t', h=F);"
	 self.init_windows! window
      else
	 f = File.open(val, "r")
	 @windows = []
	 while ln = f.gets
	    break unless /^#:/.match(ln).nil?
	    @windows << ROCWindow.new(self, ln)
	 end
	 f.close
	 @aln = Alignment.new
	 @aln.read_rocker(val)
      end
   end
   def win_at_col(col) self.windows.select{|w| (w.from<=col) and (w.to>=col)}.first end
   def in_nucl?() @nucl end
   def nucl=(nucl) @nucl=nucl end
   def refine! table
      while true
	 return false unless self.load_table! table
	 break if self._refine_iter(table)==0
      end
      return true
   end
   def _refine_iter table
      to_refine = []
      self.windows.each do |w|
	 next if w.almost_empty or w.length <= 5
	 self.rrun "acc <- w$accuracy[w$V1==#{w.from}];"
	 to_refine << w if self.rrun("ifelse(is.na(acc), 100, acc)", :float) < 95.0
      end
      n = to_refine.size
      return 0 unless n > 0
      to_refine.each do |w|
	 w1 = ROCWindow.new(self, w.from, (w.from+w.to)/2)
	 w2 = ROCWindow.new(self, (w.from+w.to)/2, w.to)
	 if w1.almost_empty or w2.almost_empty
	    n -= 1
	 else
	    @windows << w1
	    @windows << w2
	    @windows.delete w
	 end
      end
      @windows.sort!{ |x,y| x.from <=> y.from }
      n
   end
   def load_table! table, sbj=[], min_score=0
      self.rrun "x <- read.table('#{table}', sep='\\t', h=F);"
      self.rrun "x <- x[x$V1 %in% c('#{sbj.join("','")}'),];" if sbj.size > 0
      self.rrun "x <- x[x$V4 >= #{minscore.to_s},];" if min_score > 0
      Dir.mktmpdir do |dir|
         self.save(dir + "/rocker")
	 self.rrun "w <- read.table('#{dir}/rocker', sep='\\t', h=F);"
      end
      self.rrun "w <- w[!is.na(w$V5),];"
      if self.rrun("nrow(w)", :int)==0
         warn "\nWARNING: Insufficient windows with estimated thresholds.\n\n"
         return false
      end
      self.rrun <<-EOC
	 w$tp<-0; w$fp<-0; w$tn<-0; w$fn<-0;
	 for(i in 1:nrow(x)){
	    m <- x$V6[i];
	    win <- which( (m>=w$V1) & (m<=w$V2))[1];
	    if(!is.na(win)){
	       if(x$V4[i] >= w$V5[win]){
		  if(x$V5[i]==1){ w$tp[win] <- w$tp[win]+1 }else{ w$fp[win] <- w$fp[win]+1 };
	       }else{
		  if(x$V5[i]==1){ w$fn[win] <- w$fn[win]+1 }else{ w$tn[win] <- w$tn[win]+1 };
	       }
	    }
	 }
      EOC
      r.run <<-EOC
	 w$p <- w$tp + w$fp;
	 w$n <- w$tn + w$fn;
	 w$sensitivity <- 100*w$tp/(w$tp+w$fn);
	 w$specificity <- 100*w$tn/(w$fp+w$tn);
	 w$accuracy <- 100*(w$tp+w$tn)/(w$p+w$n);
	 w$precision <- 100*w$tp/(w$tp+w$fp);
      EOC
      
      return true
   end
   def init_windows!(size)
      @windows = []
      1.step(self.aln.cols,size).each { |a| @windows << ROCWindow.new(self, a, a+size-1) }
   end
   def rrun(cmd, type=nil) self.r.run cmd, type end
   def save(file)
      f = File.open(file, "w")
      f.print self.to_s
      f.close
   end
   def to_s
      o = ''
      self.windows.each{|w| o += w.to_s}
      o += self.aln.to_s
      return o
   end
end

