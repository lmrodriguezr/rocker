#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jun-08-2015
#

require 'tmpdir'

class ROCker
   #================================[ Class ]
   #@@DEFAULTS.merge!({  })
   
   #================================[ Search ]
   def search!
      raise "-k/--rocker is mandatory." if @o[:rocker].nil?
      raise "-q/--query is mandatory." if @o[:query].nil?

      # Check requirements
      puts "Testing environment." unless @o[:q]
      @o[:searchcmd] = @o[:searchcmd][@o[:search]] if @o[:searchcmd].is_a? Hash
      @o[:makedbcmd] = @o[:makedbcmd][@o[:search]] if @o[:makedbcmd].is_a? Hash
      self.bash "#{@o[:searchbins]}makeblastdb -version", "--search-bins must contain executables. Is BLAST+ installed?" if @o[:search]==:blast
      self.bash "#{@o[:searchbins]}diamond --help", "--search-bins must contain executables. Is DIAMOND installed?" if @o[:search]==:diamond

      # Run similarity search
      Dir.mktmpdir do |dir|
	 @o[:qblast] ||= "#{dir}/blast"
	 puts "Loading ROCker file: #{@o[:rocker]}." unless @o[:q]
	 data = ROCData.new @o[:rocker]
	 puts "Running similarity search." unless @o[:q]
	 puts "  * preparing database." unless @o[:q]
	 ofh = File.new("#{dir}/ref.fasta", "w")
	 ofh.print data.aln.to_seq_s
	 ofh.close
	 bash sprintf(@o[:makedbcmd], @o[:searchbins], 'prot', "#{dir}/ref.fasta", "#{dir}/ref")
	 puts "  * running similarity search." unless @o[:q]
	 bash sprintf(@o[:searchcmd], @o[:searchbins], 'blastx', @o[:query], "#{dir}/ref", @o[:qblast], @o[:thr])
	 self.filter! data
      end
   end # search!
end # ROCker

