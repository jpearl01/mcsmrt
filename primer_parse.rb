#!/usr/bin/env ruby

require 'bio'
require 'optimist'

opts = Optimist::options do
  opt :primerfile, "USEARCH results file after running the search_oligodb command", :type => :string, :short => "-p"
  opt :outfile, "Output file which will have all the information from primer matching", :type => :string, :short => "-o"
end

#TODO: add ql Query sequence length. to the usearch_local command so we don't need all the allreadsfile - DONE

#Check on inputs, error if necessary
abort("Must provide usearch local results file name") if opts[:primerfile].nil?
abort("Must provide an output file name") if opts[:outfile].nil?
abort("Usearch local results file does not exist") unless File.exist?(opts[:primerfile])

class Read
  #TODO: add default values, and define input variables to the initialize function (i.e. name="blah" in the 'new' function for the class) - DONE
  attr_accessor :name, :orientation, :fow_match, :rev_match, :note, :primers
  def initialize(name = nil, orien = nil, f = nil, r = nil, n = nil, p = [])
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
  def initialize(s = nil, e = nil, g = nil, o = nil, r = nil, l = nil)
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
#read reference_name orientation ref_star ref_end read_start read_end diffs mism gaps (TODO: check on order - DONE)

reads_hash = {}

pf = File.open(opts[:primerfile], "r")

#Initialize reads_hash with all primer matches 
pf.each do |line|
  #next unless $.<10  
  rec = line.split("\t")
  key = rec[0].split(";")[0]

  if reads_hash.has_key?(key)
    p = Primer.new(rec[5].to_i, rec[6].to_i, rec[8].to_i, rec[2], rec[1])
    reads_hash[key].primers.push(p)
  else
    reads_hash[key] = Read.new(key)
    p = Primer.new(rec[5].to_i, rec[6].to_i, rec[8].to_i, rec[2], rec[1])
    reads_hash[key].primers.push(p)
  end
end
#puts reads_hash

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

reads_hash.each do |key, value|
  if value.primers.size == 2
    if value.primers[0].reference_name != value.primers[1].reference_name
      f_ind = value.primers.find_index { |e| e.reference_name.match( /forward/ ) }
      if value.primers[f_ind].orientation == "+"
        reads_hash[key].orientation = "+"
      else
        reads_hash[key].orientation = "-"
      end
      reads_hash[key].fow_match = true
      reads_hash[key].rev_match = true
      reads_hash[key].note = "good"
    else
      reads_hash[key].orientation = "NA"
      reads_hash[key].fow_match = "NA"
      reads_hash[key].rev_match = "NA"
      reads_hash[key].note = "two_same_hits"
    end

  elsif value.primers.size > 2
    reads_hash[key].orientation = "NA"
    reads_hash[key].fow_match = "NA"
    reads_hash[key].rev_match = "NA"
    reads_hash[key].note = "more_than_2_hits"
  
  else
    reads_hash[key].orientation = "NA"
    reads_hash[key].fow_match = "NA"
    reads_hash[key].rev_match = "NA"
    reads_hash[key].note = "one_primer_hit"
  end

  if reads_hash[key].note == "good"
    f_ind = value.primers.find_index { |e| e.reference_name.match( /forward/ ) }
    r_ind = value.primers.find_index { |e| e.reference_name.match( /reverse/ ) }
    gaps_arr = [value.primers[f_ind].gaps, value.primers[r_ind].gaps]
    out_file.puts([key,
                  reads_hash[key].fow_match,
                  reads_hash[key].rev_match,
                  value.primers[f_ind].start,
                  value.primers[f_ind].end,
                  value.primers[r_ind].start,
                  value.primers[r_ind].end,
                  reads_hash[key].orientation,
                  gaps_arr.inspect,
                  reads_hash[key].note,
                  value.primers.size].join("\t")
                  )
  else
    out_file.puts([key,
                  reads_hash[key].fow_match,
                  reads_hash[key].rev_match,
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  reads_hash[key].orientation,
                  "NA",
                  reads_hash[key].note,
                  value.primers.size].join("\t")
                  )
  end

end 

#puts reads_hash



