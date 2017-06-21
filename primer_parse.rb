#!/usr/bin/env ruby

require 'bio'
require 'trollop'

opts = Trollop::options do
  opt :primerfile, "USEARCH results file after running the search_oligodb command", :type => :string, :short => "-p"
  opt :outfile, "Output file which will have all the information from primer matching", :type => :string, :short => "-o"
end

#TODO: add ql Query sequence length. to the usearch_local command so we don't need all the allreadsfile

#Check on inputs, error if necessary
abort("Must provide usearch local results file name") if opts[:primerfile].nil?
abort("Must provide an output file name") if opts[:outfile].nil?
abort("Usearch local results file does not exist") unless File.exist?(opts[:primerfile])

class Read
  #TODO: add default values, and define input variables to the initialize function (i.e. name="blah" in the 'new' function for the class)
  attr_accessor :name, :orientation, :fow_match, :rev_match, :note, :primers
  def initialize(name, orien, f, r, n, p)
    @name = name
    @orientation = orien
    @fow_match = f
    @rev_match = r
    @note = n
    @primers = p
  end
end

class Primer
  attr_accessor :start, :end, :gaps, :orientation, :reference_name
  def initialize(s, e, g, o, r)
    #it would also be nice to check and make sure each of the things being entered is checked to make sure it is the type of value we are expecting
    @start = s
    @end = e
    @gaps = g
    @orientation = o
    @reference_name = r
  end
end


#Parse primer match file
#File Structure:
#read refernce_name orientation ref_star ref_end read_start read_end diffs mism gaps (TODO: check on order)

reads_hash = {}

pf = File.open(opts[:primerfile], "r")

pf.each do |line|
  next unless $.<10  
  rec = line.split("\t")

  if reads_hash.has_key?(rec[0])
    
  else
    reads_hash[rec[0]] = Read.new(rec[0], nil, nil, nil, nil, [])
    #TODO check on the order in the primer_map.txt file
    p = Primer.new(rec[5], rec[6], rec[10], rec[2], rec[1])
    reads_hash[rec[0]].primers.push(p)
  end


end

puts reads_hash

out_file = File.open(opts[:outfile], "w")
out_file.puts(["read_name", 
               "forward_primer_match",
               "reverse_primer_match",
               "forward_primer_start",
               "forward_primer_end", 
               "reverse_primer_start", 
               "reverse_primer_end", 
               "read_orientation",
               "gaps", 
               "primer_note", 
               "num_of_primer_hits"].join("\t")
              )
