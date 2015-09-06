#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Sep-06-2015
#

class BlastHit
   attr_reader :sbj, :sfrom, :sto, :bits, :istrue, :isfalse, :midpoint
   # Initialize from BLAST using new(ln,aln),
   # initialize from TABLE using new(ln)
   def initialize(ln, aln=nil)
      l = ln.chomp.split(/\t/)
      if aln.nil?
	 @sbj	= l[0]
	 @sfrom	= l[1].to_i
	 @sto	= l[2].to_i
	 @bits	= l[3].to_f
	 @istrue = l[4]=='1'
	 @istrue = l[4]=='-1'
	 @midpoint = l[5].to_i
      else
	 s = aln.seq(l[1])
	 return nil if s.nil?
	 @sbj	= s.id
	 a	= s.pos2col(l[8].to_i)
	 b 	= s.pos2col(l[9].to_i)
	 @sfrom	= [a,b].min
	 @sto	= [a,b].max
	 @bits	= l[11].to_f
	 @istrue = ! /@%/.match(l[0]).nil?
	 @isfalse = ! /@\$/.match(l[0]).nil?
	 @midpoint = s.pos2col(((l[8].to_f+l[9].to_f)/2).ceil)
      end
   end
   def to_s
      self.sbj.nil? ? "" :
	 [sbj, sfrom.to_s, sto.to_s, bits.to_s,
	    istrue ? "1" : (isfalse ? "-1" : "0"), midpoint].join("\t") + "\n"
   end
end

