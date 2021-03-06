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
    readlengths = {}
    exp_readlen = 0
    unless @o[:lencorr].nil?
      @o[:lencorr_max] ||= 0.4
      raise "Unsigned length in model, please re-compile model to use -L" if
        data.signatures[:l].nil?
      exp_readlen = data.signatures[:l].to_i
      File.open(@o[:lencorr], 'r') do |fh|
        k = nil
        fh.each_line do |ln|
          if ln =~ /^>(\S+)/
            k = $1
            readlengths[k] = 0
          else
            readlengths[k] += ln.chomp.size
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
        next if bh.sbj.nil? # <- When the hit is not against a known target
        bs = @o[:lencorr].nil? ? bh.bits :
          correct_bs(bh, readlengths[bh.qry], exp_readlen,
            @o[:lencorr_max], @o[:lencorr_pen])
        oh.print ln if not(bh.sfrom.nil?) and
          bs >= data.win_at_col(bh.midpoint).thr
      end
    end
    oh.close
  end # filter!

  def correct_bs(bh, readlen, exp_readlen, max_corr, penalty)
    bs = bh.bits
    return bs if @o[:lencorr].nil? or readlen.nil? or readlen >= exp_readlen
    bits_per_aa = bs.to_f / readlen
    if penalty.nil? or penalty > 1.0
      extra = [0.0, readlen * (max_corr + 1.0) - exp_readlen].max
      max_tri = max_corr * readlen * bits_per_aa / 2
      tanTheta = max_corr > 0.0 ? bits_per_aa / (max_corr * readlen) : 0.0
      extra_tri = extra * extra * tanTheta / 2
      return bs + (max_tri - extra_tri)
    else
      miss = [exp_readlen - readlen, max_corr * readlen].min
      return bs + (bits_per_aa * miss * (1.0 - penalty))
    end
  end
end # ROCker

