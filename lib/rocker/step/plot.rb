#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Sep-09-2015
#

class ROCker
   #================================[ Class ]
   @@DEFAULTS.merge!({
      color:false, gformat:"pdf", width:9, height:9, impact:false,
      transparency:true, sbj:[], tag_negatives:false
   })
   
   #================================[ Search ]
   def plot!
      raise "-k/--rocker is mandatory." if o[:rocker].nil?
      if @o[:table].nil?
	 raise "-t/--table is mandatory unless -b is provided." if
	    @o[:blast].nil?
	 @o[:table] = "#{@o[:blast]}.table"
      end
      raise "-b/--blast is mandatory unless -t exists." if
	 @o[:blast].nil? and not File.exist? @o[:table]

      puts "Testing environment." unless @o[:q]
      bash "echo '' | #{@o[:r]} --vanilla", "-r/--path-to-r must be " +
	 "executable. Is R installed?"

      # Source files
      puts "Reading files." unless @o[:q]
      puts "  * loding ROCker file: #{@o[:rocker]}." unless @o[:q]
      data = ROCData.new @o[:rocker]
      if File.exist? @o[:table]
	 puts "  * reusing existing file: #{@o[:table]}." unless @o[:q]
      else
	 puts "  * generating table: #{@o[:table]}." unless @o[:q]
	 blast2table(@o[:blast], @o[:table], data.aln, @o[:minscore])
      end

      # Matches (middle panel)
      puts "Plotting matches." unless @o[:q]
      extra = @o[:gformat]=="pdf" ? "" : ", units='in', res=300"
      @o[:gout] ||= "#{@o[:rocker]}.#{@o[:gformat]}"
      data.rrun "#{@o[:gformat]}('#{@o[:gout]}', #{@o[:width]}, " +
	 "#{@o[:height]}#{extra});"
      data.rrun "layout(c(2,1,3), heights=c(2-1/#{data.aln.size},3,1));"
      some_thr = data.load_table! @o[:table], @o[:sbj], @o[:minscore]
      data.rrun "par(mar=c(0,4,0,0.5)+.1);"
      data.rrun "plot(1, t='n', xlim=c(0.5,#{data.aln.cols}+0.5), " +
	 "ylim=range(x$V4)+c(-0.04,0.04)*diff(range(x$V4)), xlab='', " +
	 "ylab='Bit score', xaxs='i', xaxt='n');"
      data.rrun "noise <- runif(ncol(x),-.2,.2)"
      data.rrun "hit.col <- ifelse(x$V5==1, " +
	 "rgb(0,0,.5,#{@o[:transparency] ? ".2" : "1"}), " +
	 "rgb(.5,0,0,#{@o[:transparency] ? ".2" : "1"}))"
      data.rrun "hit.col[ x$V5==-1 ] <- " +
	 "rgb(0.722,0.722,0,#{@o[:transparency] ? ".2" : "1"})" if
	 @o[:tag_negatives]
      data.rrun "arrows(x0=x$V2, x1=x$V3, y0=x$V4+noise, lty=1, col=hit.col, " +
	 "length=0);"
      data.rrun "points(x$V6, x$V4+noise, col=hit.col, pch=19, cex=1/4);"
      
      # Windows (middle panel)
      puts "Plotting windows." unless @o[:q]
      if some_thr
	 data.rrun "arrows(x0=w$V1, x1=w$V2, y0=w$V5, lwd=2, length=0)"
	 data.rrun "arrows(x0=w$V2[-nrow(w)], x1=w$V1[-1], " +
	    "y0=w$V5[-nrow(w)], y1=w$V5[-1], lwd=2, length=0)"
      end
      data.rrun "legend('bottomright', legend=c('Match span'," +
	 "'Match mid-point','Reference (+)'," +
	 "#{"'Reference (-)'," if @o[:tag_negatives]}'Non-reference'), " +
	 "lwd=c(1,NA,1,1,1), pch=c(NA,19,19,19,19), ncol=5, bty='n', " +
	 "col=c('black','black','darkblue'," +
	 "#{"rgb(.722,.722,0)," if @o[:tag_negatives]}'darkred'))"

      # Alignment (top panel)
      puts "Plotting alignment." unless @o[:q]
      data.rrun "par(mar=c(0,4,0.5,0.5)+0.1);"
      data.rrun "plot(1, t='n', xlim=c(0,#{data.aln.cols}), " +
	 "ylim=c(1,#{data.aln.seqs.size}), xlab='', ylab='Alignment', " +
	 "xaxs='i', xaxt='n', yaxs='i', yaxt='n', bty='n');"
      i = 0
      data.rrun "clr <- rainbow(26, v=1/2, s=3/4);" if @o[:color]
      data.aln.seqs.values.each do |s|
         color = (s.aln.split(//).map do |c|
	    c=="-" ? "'grey80'" :
	       (@o[:sbj].include?(s.id) ? "'red'" :
	       (@o[:color] ? "clr[#{c.ord-64}]" :
	       "'black'"))
	 end.join(","))
	 data.rrun "rect((1:#{data.aln.cols-1})-0.5, " +
	    "rep(#{i}, #{data.aln.cols-1}), (1:#{data.aln.cols-1})+0.5, " +
	    "rep(#{i+1}, #{data.aln.cols-1}), col=c(#{color}), border=NA);"
	 i += 1
      end

      # Statistics (bottom panel)
      puts "Plotting statistics." unless @o[:q]
      data.rrun "par(mar=c(5,4,0,0.5)+.1);"
      data.rrun "plot(1, t='n', xlim=c(0,#{data.aln.cols}), " +
	 "ylim=c(#{@o[:ylim].nil? ? (@o[:impact] ? "-2,.1" : "50,100") :
	    @o[:ylim]}), xlab='Alignment position (amino acids)', " +
	 "ylab='Precision',xaxs='i');"
      if some_thr
	 sn = data.rrun "100*sum(w$tp)/(sum(w$tp)+sum(w$fn))", :float
	 sp = data.rrun "100*sum(w$tn)/(sum(w$fp)+sum(w$tn))", :float
	 ac = data.rrun "100*(sum(w$tp)+sum(w$tn))/(sum(w$p)+sum(w$n))", :float
	 unless @o[:q]
	    puts "  * sensitivity: #{sn}%"
	    puts "  * specificity: #{sp}%"
	    puts "  * accuracy: #{ac}%"
	 end
	 data.rrun "pos <- (w$V1+w$V2)/2"
	 if @o[:impact]
	    data.rrun "lines(pos[!is.na(w$specificity)], " +
	       "(w$specificity[!is.na(w$specificity)]-#{sp})*" +
		  "w$tp[!is.na(w$specificity)]/sum(w$tp), " +
	       "col='darkred', lwd=2, t='o', cex=1/3, pch=19);"
	    data.rrun "lines(pos[!is.na(w$sensitivity)], " +
	       "(w$sensitivity[!is.na(w$sensitivity)]-#{sn})*" +
		  "w$tn[!is.na(w$sensitivity)]/sum(w$tn), " +
	       "col='darkgreen', lwd=2, t='o', cex=1/3, pch=19);"
	    data.rrun "lines(pos[!is.na(w$accuracy)], " +
	       "(w$accuracy[!is.na(w$accuracy)]-#{ac})*" +
		  "(w$tp+w$tn)[!is.na(w$accuracy)]/sum(c(w$tp, w$tn)), " +
	       "col='darkblue', lwd=2, t='o', cex=1/3, pch=19);"
	 else
	    data.rrun "lines(pos[!is.na(w$specificity)], " +
	       "w$specificity[!is.na(w$specificity)], col='darkred', " +
	       "lwd=2, t='o', cex=1/3, pch=19);"
	    data.rrun "lines(pos[!is.na(w$sensitivity)], " +
	       "w$sensitivity[!is.na(w$sensitivity)], col='darkgreen', " +
	       "lwd=2, t='o', cex=1/3, pch=19);"
	    data.rrun "lines(pos[!is.na(w$accuracy)], " +
	       "w$accuracy[!is.na(w$accuracy)], col='darkblue', lwd=2, " +
	       "t='o', cex=1/3, pch=19);"
	 end
      end
      data.rrun "legend('bottomright', " +
	 "legend=c('Specificity','Sensitivity','Accuracy'), lwd=2, " +
	 "col=c('darkred','darkgreen','darkblue'), ncol=3, bty='n')"
      data.rrun "dev.off();"
   end # plot!
end # ROCker

