require 'trollop'

#ruby retrieve_singletons.rb -f reads_of_insert_without_barcodes_singletons.fasta -p reads_of_insert_without_barcodes_singletons_primer_map.txt -w reads_of_insert_without_barcodes_singletons_primer_map_forward.fasta -r reads_of_insert_without_barcodes_singletons_primer_map_reverse.fasta -o retrieved_singleton_seqs.fasta

opts = Trollop::options do
	opt :primerfile, "USEARCH results file after running the search_oligodb command", :type => :string, :short => "-p"
	opt :outfile, "Output file which will have all the information from primer matching", :type => :string, :short => "-o"
end
#puts opts

abort("!!!!The primer matching results file does not exist!!!!") if !File.exists?(opts[:primerfile])

primer_file = File.open(opts[:primerfile], "r")
out_file = File.open(opts[:outfile], "w")
out_file.puts("read_name\tforward_primer_match\treverse_primer_match\tforward_primer_start\tforward_primer_end\treverse_primer_start\treverse_primer_end\tread_orientation\tprimer_note")

class Primer_info
  attr_accessor :f_primer_matches, :r_primer_matches, :f_primer_start, :f_primer_end, :r_primer_start, :r_primer_end, :read_orientation, :primer_note

  def initialize
    @f_primer_matches = false
    @r_primer_matches = false
    @f_primer_start = 0
    @f_primer_end = 0
    @r_primer_start = 0
    @r_primer_end = 0
    @read_orientation = ""
    @primer_note = ""
  end
end


def parse_primer_file (primer_file)
	# Storing the record header as the key and primer strand and sequence strand as the value. 
	record_hash = {}
    primer_file.each do |line|
    	line_array = line.split("\t")
      	key = line_array[0].split(";")[0]
      	if record_hash.has_key?(key)
        	record_hash[key].push([line_array[1], line_array[2], line_array[6].chomp.to_i, line_array[7].chomp.to_i])
	    else
	        record_hash[key] = [[line_array[1], line_array[2], line_array[6].chomp.to_i, line_array[7].chomp.to_i]]
	    end
	end
	#puts record_hash 
	return record_hash                                                                                                                                     
end 

def print_info (record_hash, primer_hash, out_file)
	record_hash.each do |key, value|
		primer_hash[key] = Primer_info.new
  		# When size of value == 2, exactly one primer pair is matched
  		if value.size == 2
  			if value[0][0].eql?("forward") and value[1][0].eql?("forward")
  				primer_hash[key].f_primer_matches = true
  				primer_hash[key].r_primer_matches = false
  				primer_hash[key].f_primer_start = "NA"
  				primer_hash[key].f_primer_end = "NA"
  				primer_hash[key].r_primer_start = "NA"
  				primer_hash[key].r_primer_end = "NA"
  				primer_hash[key].read_orientation = "NA"
  				primer_hash[key].primer_note = "two_forward_primer_hits" 
  		
  			elsif value[0][0].eql?("reverse") and value[1][0].eql?("reverse")
  				primer_hash[key].f_primer_matches = false
  				primer_hash[key].r_primer_matches = true
  				primer_hash[key].f_primer_start = "NA"
  				primer_hash[key].f_primer_end = "NA"
  				primer_hash[key].r_primer_start = "NA"
  				primer_hash[key].r_primer_end = "NA"
  				primer_hash[key].read_orientation = "NA"
  				primer_hash[key].primer_note = "two_reverse_primer_hits"
  		
  			elsif value[0][0].eql?("forward") and value[1][0].eql?("reverse")
  				primer_hash[key].f_primer_matches = true
  				primer_hash[key].r_primer_matches = true
  				primer_hash[key].f_primer_start = value[0][2].to_i
  				primer_hash[key].f_primer_end = value[0][3].to_i
  				primer_hash[key].r_primer_start = value[1][2].to_i
  				primer_hash[key].r_primer_end = value[1][3].to_i
          primer_hash[key].primer_note = "NA" 
  				if value[0][1] == "+"
  					primer_hash[key].read_orientation = "+"
  				else
  					primer_hash[key].read_orientation = "-"
  				end

  			elsif value[1][0].eql?("forward") and value[0][0].eql?("reverse")  
  				primer_hash[key].f_primer_matches = true
  				primer_hash[key].r_primer_matches = true
  				primer_hash[key].f_primer_start = value[1][2].to_i
  				primer_hash[key].f_primer_end = value[1][3].to_i
  				primer_hash[key].r_primer_start = value[0][2].to_i
  				primer_hash[key].r_primer_end = value[0][3].to_i
          primer_hash[key].primer_note = "NA" 
  				if value[1][1] == "+"
  					primer_hash[key].read_orientation = "+"
  				else
  					primer_hash[key].read_orientation = "-"
  				end
  			end

      elsif value.size > 2
        if value[0][0].eql?("reverse") and value[1][0].eql?("forward")
          primer_hash[key].r_primer_matches = true 
          primer_hash[key].f_primer_matches = true 
          primer_hash[key].f_primer_start   = value[1][2].to_i
          primer_hash[key].f_primer_end     = value[1][3].to_i
          primer_hash[key].r_primer_start   = value[0][2].to_i
          primer_hash[key].r_primer_end     = value[0][3].to_i
          if value[1][1] == "+"
            primer_hash[key].read_orientation = "+"
          else
            primer_hash[key].read_orientation = "-"
          end
          primer_hash[key].primer_note = "multiple_primers_hits"

        elsif value[0][0].eql?("forward") and value[1][0].eql?("reverse")
          primer_hash[key].r_primer_matches = true 
          primer_hash[key].f_primer_matches = true 
          primer_hash[key].f_primer_start   = value[0][2].to_i
          primer_hash[key].f_primer_end     = value[0][3].to_i
          primer_hash[key].r_primer_start   = value[1][2].to_i
          primer_hash[key].r_primer_end     = value[1][3].to_i
          if value[0][1] == "+"
            primer_hash[key].read_orientation = "+"
          else
            primer_hash[key].read_orientation = "-"
          end
          primer_hash[key].primer_note = "multiple_primers_hits"

        else
          value.each do |rec|
            primer_hash[key].r_primer_matches = true if rec[0].eql?("reverse")
            primer_hash[key].f_primer_matches = true if rec[0].eql?("forward")
          end
          primer_hash[key].f_primer_start   = "NA"
          primer_hash[key].f_primer_end     = "NA"
          primer_hash[key].r_primer_start   = "NA"
          primer_hash[key].r_primer_end     = "NA"
          primer_hash[key].read_orientation = "NA"
          primer_hash[key].primer_note = "multiple_primers_hits"
        end

  		# When size of value == 1, its a singleton
  		elsif value.size == 1
  			if value[0][0].eql?("forward")
  				primer_hash[key].f_primer_matches = true
  				primer_hash[key].f_primer_start = value[0][2].to_i
  				primer_hash[key].f_primer_end = value[0][3].to_i
  				if value[0][1] == "+"
  					primer_hash[key].read_orientation = "+"
  				else
  					primer_hash[key].read_orientation = "-"
  				end
  				primer_hash[key].primer_note = "reverse_missing_singleton"
  			end

  			if value[0][0].eql?("reverse")
  				primer_hash[key].r_primer_matches = true
  				primer_hash[key].r_primer_start = value[0][2].to_i
  				primer_hash[key].r_primer_end = value[0][3].to_i
  				if value[0][1] == "+"
  					primer_hash[key].read_orientation = "-"
  				else
  					primer_hash[key].read_orientation = "+"
  				end
  				primer_hash[key].primer_note = "forward_missing_singleton"
  			end
  		end
  		
  		out_file.puts("#{key}\t#{primer_hash[key].f_primer_matches}\t#{primer_hash[key].r_primer_matches}\t#{primer_hash[key].f_primer_start}\t#{primer_hash[key].f_primer_end}\t#{primer_hash[key].r_primer_start}\t#{primer_hash[key].r_primer_end}\t#{primer_hash[key].read_orientation}\t#{primer_hash[key].primer_note}")

  	end
  	#puts primer_hash["m151002_181152_42168_c100863312550000001823190302121650_s1_p0/85/ccs"].inspect #test for singleton forward
	#puts all_reads_hash["m151002_181152_42168_c100863312550000001823190302121650_s1_p0/66/ccs"].inspect #test for more than 2 primer hits
	#puts primer_hash.length
end

primer_hash = {}
record_hash = parse_primer_file(primer_file)
#puts record_hash.length
print_info(record_hash, primer_hash, out_file)