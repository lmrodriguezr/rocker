#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Sep-07-2015
#

class ROCker
   #================================[ Class ]
   @@DEFAULTS.merge!({:refine=>true, :win=>20, :minscore=>0})
   
   #================================[ Compile ]
   def compile!
      raise "-a/--alignment is mandatory." if @o[:aln].nil?
      raise "-a/--alignment must exist." unless File.exist? @o[:aln]
      raise "-l/--readlen is mandatory." if @o[:readlen].nil?
      if @o[:table].nil?
	 raise "-T/--table is mandatory unless -b is provided." if
	    @o[:blast].nil? or not File.exist? @o[:blast]
	 @o[:table] = "#{@o[:blast]}.table"
      else
	 @o[:reuse] = true
      end
      raise "-b/--blast is mandatory unless -t exists." if
	 @o[:blast].nil? and not File.exist? @o[:table]
      raise "-k/--rocker is mandatory." if @o[:rocker].nil?

      puts "Testing environment." unless @o[:q]
      bash("echo '' | #{@o[:r]} --vanilla",
	 "-r/--path-to-r must be executable. Is R installed?")
      bash("echo \"library('pROC')\" | #{@o[:r]} --vanilla",
	 "Please install the 'pROC' library for R first.")

      puts "Reading files." unless @o[:q]
      puts "  * loading alignment: #{@o[:aln]}." unless @o[:q]
      aln = Alignment.new
      aln.read_fasta @o[:aln]
      
      if @o[:reuse] and File.exist? @o[:table]
	 puts "  * reusing existing file: #{@o[:table]}." unless @o[:q]
      else
	 puts "  * generating table: #{@o[:table]}." unless @o[:q]
	 blast2table(@o[:blast], @o[:table], aln, @o[:minscore])
      end

      puts "Analyzing data." unless @o[:q]
      puts "  * computing windows." unless @o[:q]
      data = ROCData.new(@o[:table], aln, @o[:win])
      data.nucl = @o[:nucl]
      if @o[:refine]
	 puts "  * refining windows." unless @o[:q]
	 warn "Insufficient hits to refine results." unless
	    data.refine! @o[:table]
      end
      puts "  * saving ROCker file: #{@o[:rocker]}." unless @o[:q]
      data.save(@o[:rocker], l: @o[:readlen])
   end # compile!
end # ROCker

