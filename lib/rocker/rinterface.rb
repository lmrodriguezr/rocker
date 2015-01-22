#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jan-22-2015
#

class RInterface
   @@R_BIN = "R"
   def RInterface.R_BIN=(rbin) @@R_BIN=rbin end
   attr_reader :handler
   def initialize
      @handler = IO.popen("#{@@R_BIN} --slave 2>&1", "w+")
   end
   def run(cmd, type=nil)
      @handler.puts cmd
      @handler.puts "cat('---FIN---\n')"
      o = ""
      while true
         l = @handler.gets
	 raise "R failed on command:\n#{cmd}\n\nError:\n#{o}" if l.nil?
	 break unless /^---FIN---/.match(l).nil?
	 o += l
      end
      o.chomp!
      case type
      when :float
	 /^\s*\[1\]\s+([0-9\.Ee+-]+|Inf).*/.match(o).nil? and raise "R error: expecting float, got #{o}"
	 return Float::INFINITY if $1=='Inf'
	 return $1.to_f
      when :int
	 /^\s*\[1\]\s+([0-9\.Ee+-]+).*/.match(o).nil? and raise "R error: expecting integer, got #{o}"
	 return $1.to_i
      else
	 return o
      end
   end
end

