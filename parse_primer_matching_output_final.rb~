require 'bio'
require 'trollop'

opts = Trollop::options do
	opt :fastqfile, "Fastq file with all the records", :type => :string, :short => "-f"
	opt :primerfile, "Primer matched file which is the output of the oligo_db command", :type => :string, :short => "-p"
	opt :outfile, "Output file to write the oriented sequences in", :type => :string, :short => "-o"
	opt :outfile_2, "Output file to write the oriented singleton sequences in", :type => :string, :short => "-s"
	opt :outfile_3, "Output file with the primer map details for singletons", :type => :string, :short => "-m"
end
#puts opts
abort("!!!!The primer matched output file (required for parsing the primer matching output) does not exist!!!!") if !File.exists?(opts[:primerfile])
abort("!!!!The fastq file (with all the sequences which were primer matched) does not exist!!!!") if !File.exists?(opts[:fastqfile])

primer_file = File.open(opts[:primerfile], "r")
fastq_file = Bio::FlatFile.auto(opts[:fastqfile])
out_file = File.open(opts[:outfile], "w")
out_file_2 = File.open(opts[:outfile_2], "w")
out_file_3 = File.open(opts[:outfile_3], "w")

record_hash = {}
primer_file.each do |line|
	line_array = line.split("\t")
	if record_hash.has_key?(line_array[0])
		record_hash[line_array[0]].push([line_array[1], line_array[2], line_array[6].chomp.to_i, line_array[7].chomp.to_i])	
  	else
    	record_hash[line_array[0]] = [[line_array[1], line_array[2], line_array[6].chomp.to_i, line_array[7].chomp.to_i]]
 	end
end
#puts record_hash

singletons = {}
to_complement_hash = {}
record_hash.each do |key, value|
	if value.size == 2
    		if value[0][0].eql?("reverse") && value[0][1].eql?("-") && value[1][0].eql?("forward") && value[1][1].eql?("+")
      			to_complement_hash[key] = "Oriented" 
    		elsif value[0][0].eql?("reverse") && value[0][1].eql?("+") && value[1][0].eql?("forward") && value[1][1].eql?("-")
      			to_complement_hash[key] = "Orient"
    		elsif value[0][0].eql?("forward") && value[0][1].eql?("+") && value[1][0].eql?("reverse") && value[1][1].eql?("-")
      			to_complement_hash[key] = "Oriented"
    		elsif value[0][0].eql?("forward") && value[0][1].eql?("-") && value[1][0].eql?("reverse") && value[1][1].eql?("+")
      			to_complement_hash[key] = "Orient"
    		end	 
  	elsif value.size == 1
		out_file_3.puts(key+"\t"+value[0][0]+"\t"+value[0][1]+"\t"+value[0][2].to_s+"\t"+value[0][3].to_s)
    		if value[0][0].eql?("reverse") && value[0][1].eql?("-") 
      			singletons[key] = "Oriented"
    		elsif value[0][0].eql?("reverse") && value[0][1].eql?("+") 
      			singletons[key] = "Orient"
    		elsif value[0][0].eql?("forward") && value[0][1].eql?("+") 
      			singletons[key] = "Oriented"
    		elsif value[0][0].eql?("forward") && value[0][1].eql?("-") 
      			singletons[key] = "Orient"
   		end	
  	end
end

fastq_file.each do |entry|
	definition = entry.definition.split(" ")[0]
	if to_complement_hash[definition].eql?("Oriented")
		if record_hash[definition][0][0].eql?("reverse")
			#puts "ONE #{record_hash[definition][1][3]} to #{record_hash[definition][0][2]}"
			out_file.puts("@"+entry.definition.split(" ")[0])
			out_file.puts(entry.naseq.upcase[record_hash[definition][1][3]..record_hash[definition][0][2]])
    			out_file.puts("+")
    			out_file.puts(entry.quality_string[record_hash[definition][1][3]..record_hash[definition][0][2]])
		else
			#puts "TWO #{record_hash[definition][0][3]} to #{record_hash[definition][1][2]}"
			out_file.puts("@"+entry.definition.split(" ")[0])
			out_file.puts(entry.naseq.upcase[record_hash[definition][0][3]..record_hash[definition][1][2]])
    			out_file.puts("+")
    			out_file.puts(entry.quality_string[record_hash[definition][0][3]..record_hash[definition][1][2]])
		end
  	elsif to_complement_hash[definition].eql?("Orient")
		if record_hash[definition][0][0].eql?("reverse")
			#puts "ONE #{record_hash[definition][0][3]} to #{record_hash[definition][1][2]}"
    			out_file.puts("@"+entry.definition.split(" ")[0])
    			out_file.puts(entry.naseq.complement.upcase[record_hash[definition][0][3]..record_hash[definition][1][2]])
    			out_file.puts("+")
    			out_file.puts(entry.quality_string[record_hash[definition][0][3]..record_hash[definition][1][2]])
		else
			#puts "TWO #{record_hash[definition][1][3]} to #{record_hash[definition][0][2]}"
			out_file.puts("@"+entry.definition.split(" ")[0])
    			out_file.puts(entry.naseq.complement.upcase[record_hash[definition][1][3]..record_hash[definition][0][2]])
    			out_file.puts("+")
    			out_file.puts(entry.quality_string[record_hash[definition][1][3]..record_hash[definition][0][2]])
		end
	end
  	if singletons[definition].eql?("Orient")
    		out_file_2.puts("@"+entry.definition.split(" ")[0])
    		out_file_2.puts(entry.naseq.complement.upcase)
		out_file_2.puts("+")
    		out_file_2.puts(entry.quality_string)
  	elsif singletons[definition].eql?("Oriented")
    		out_file_2.puts("@"+entry.definition.split(" ")[0])
    		out_file_2.puts(entry.naseq.upcase)
		out_file_2.puts("+")
    		out_file_2.puts(entry.quality_string)		
  	end 	
end
#puts no_category_hash_2, no_category_hash_2.length


