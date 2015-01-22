#
# @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
# @author Luis (Coto) Orellana
# @license artistic license 2.0
# @update Jan-22-2015
#

require 'rocker/sequence'


class Alignment
   attr_reader :seqs, :cols
   def initialize
      @seqs = {}
   end
   def read_fasta(file) self.read_file(file, false) end
   def read_rocker(file) self.read_file(file, true) end
   def read_file(file, is_rocker)
      f = File.open(file, 'r')
      id = nil
      sq = ""
      while ln = f.gets
	 if is_rocker
	    next if /^#:(.*)/.match(ln).nil?
	    ln = $1
	 end
	 m = /^>(\S+)/.match(ln)
	 if m.nil?
	    sq += ln
	 else
	    self << Sequence.new(id, sq) unless id.nil?
	    id = m[1]
	    sq = ""
	 end
      end
      self << Sequence.new(id, sq) unless id.nil?
   end
   def <<(seq)
      @seqs[seq.id] = seq
      @cols = seq.cols if self.cols.nil?
      raise "Aligned sequence #{seq.id} has a different length (#{seq.cols} vs #{self.cols})" unless seq.cols == self.cols
   end
   def get_gis
      regexps = [/^gi\|(\d+)\|/, /^(\d+)\|/, /^(\d+)$/, /^gi\|(\d+)$/, /\|gi\|(\d+)\|/, /\|gi\|(\d+)$/]
      gis = []
      self.seqs.keys.each do |id|
	 gi = nil
	 regexps.each do |regexp|
	    unless regexp.match(id).nil?
	       gi = $1
	       break
	    end
	 end
	 gis << gi unless gi.nil?
      end
      gis
   end
   def seq(id) @seqs[id] end
   def size() self.seqs.size end
   def to_seq_s() self.seqs.values.map{|s| s.to_seq_s}.join + "\n" end
   def to_s() self.seqs.values.map{|s| s.to_s}.join + "\n" end
end

