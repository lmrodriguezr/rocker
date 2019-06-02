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
    corr = {}
    readlen = 0
    unless @o[:lencorr].nil?
      raise "Unsigned length in model, please re-compile model to use -L" if
        data.signatures[:l].nil?
      readlen = data.signatures[:l].to_i
      File.open(@o[:lencorr], 'r') do |fh|
        k = nil
        fh.each_line do |ln|
          if ln =~ /^>(\S+)/
            k = $1
            corr[k] = 0
          else
            corr[k] += ln.chomp.size
          end
        end
      end
    end

    # Filter similarity search
    puts "Filtering similarity search: #{@o[:qblast]}." unless @o[:q]
    oh = File.open(@o[:oblast], 'w')
    File.open(@o[:qblast], 'r') do |ih|
      ih.each_line do |ln|
        bh = BlastHit.new(ln, data.aln)
        bs = bh.bits
        unless @o[:lencorr].nil?
          corrlen = [corr[bh.qry].to_i, 0.6 * readlen].max
          bs = bs * readlen / corrlen if corrlen < readlen
        end
        oh.print ln if
          not(bh.sfrom.nil?) and bs >= data.win_at_col(bh.midpoint).thr
      end
    end
    oh.close
  end # filter!
end # ROCker

