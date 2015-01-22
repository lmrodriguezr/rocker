#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jan-22-2015
#

class Sequence
   attr_reader :id, :seq, :aln
   def initialize(id, aln)
      @id = id
      @aln = aln.gsub(/[-\.]/,'-').gsub(/[^A-Za-z-]/, '').upcase
      @seq = aln.gsub(/[^A-Za-z]/, '').upcase
   end
   def pos2col(pos)
      col = 0
      self.aln.split(//).each do |c|
	 col+=1
	 pos-=1 unless c=='-'
	 return col if pos==0
      end
      col
   end
   def col2pos(col)
      pos = 1
      self.aln.split(//).each do |c|
         col-=1
	 pos+=1 unless c=='-'
	 return pos if col==0
      end
      pos
   end
   def cols() self.aln.length end
   def length() self.seq.length end
   def to_seq_s() ">#{self.id}\n#{self.seq}\n" end
   def to_s() "#:>#{self.id}\n#:#{self.aln}\n" end
end

