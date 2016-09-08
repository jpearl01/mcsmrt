require 'bio'
require 'trollop'

#ruby retrieve_singletons.rb -f reads_of_insert_without_barcodes_singletons.fasta -p reads_of_insert_without_barcodes_singletons_primer_map.txt -w reads_of_insert_without_barcodes_singletons_primer_map_forward.fasta -r reads_of_insert_without_barcodes_singletons_primer_map_reverse.fasta -o retrieved_singleton_seqs.fasta

opts = Trollop::options do
	opt :fastqfile, "Fastq file with all the singletons", :type => :string, :short => "-f"
	opt :primerfile, "Primer matched file which is the output of the oligo_db command for singletons", :type => :string, :short => "-p"
	opt :outfile_fow, "Fastq file with sequences which matched to the full forward primer", :type => :string, :short => "-w"
	opt :outfile_rev, "Fastq file with sequences which matched to the full reverse primer", :type => :string, :short => "-r"
	opt :outfile_fow_retrieved, "Fastq file with sequences which were retrieved from full forward primer match", :type => :string, :short => "-x"
	opt :outfile_rev_retrieved, "Fastq file with sequences which were retrieved from full reverse primer match", :type => :string, :short => "-s"
end
#puts opts

abort("!!!!The primer matched output file (required for parsing the primer matching output for singletons) does not exist!!!!") if !File.exists?(opts[:primerfile])
abort("!!!!The fastq file (with all the sequences which were primer matched to half the fow or rev primer) does not exist!!!!") if !File.exists?(opts[:fastqfile])

primer_file = File.open(opts[:primerfile], "r")
fastq_file = Bio::FlatFile.auto(opts[:fastqfile])
out_file_fow = File.open(opts[:outfile_fow], "w")
out_file_rev = File.open(opts[:outfile_rev], "w")
out_file_fow_retreived = File.open(opts[:outfile_fow_retrieved], "w")
out_file_rev_retreived = File.open(opts[:outfile_rev_retrieved], "w")

# Get the primers mapping to forward and reverse separately in a hash
forward_hash = {}
reverse_hash = {}
primer_file.each do |line|
	line_split = line.split("\t")
	key_1 = line_split[0]
  	value_1 = line_split[1]
  	if value_1.eql?("forward")
    	forward_hash[key_1] = value_1
  	else
    	reverse_hash[key_1] = value_1
  	end
end

# Create different fasta files for the seqs matching to the 2 primers
fastq_hash = {}
fastq_file.each do |entry|
	definition = entry.definition.split(" ")[0]
	fastq_hash[entry.definition] = [entry.naseq.upcase, entry.quality_string]
	if forward_hash[definition].eql?("forward")
		out_file_fow.puts("@"+entry.definition.split(" ")[0])
		out_file_fow.puts(entry.naseq.upcase)
		out_file_fow.puts("+")
		out_file_fow.puts(entry.quality_string)
	end
	if reverse_hash[definition].eql?("reverse")
		out_file_rev.puts("@"+entry.definition.split(" ")[0])
		out_file_rev.puts(entry.naseq.upcase)
		out_file_rev.puts("+")
		out_file_rev.puts(entry.quality_string)
	end
end
out_file_fow.close
out_file_rev.close

# Run usearch on these 2 fasta files with half the primer which was not getting matched
fow_basename = File.basename(opts[:outfile_fow], ".fastq")
rev_basename = File.basename(opts[:outfile_rev], ".fastq")
if !File.zero?("#{fow_basename}.fastq")
	`usearch -search_oligodb #{fow_basename}.fastq -db primer_half_rev.fasta -strand both -userout #{fow_basename}_map2.txt -userfields query+target+qstrand+diffs+tlo+thi+qlo+qhi`
end
if !File.zero?("#{rev_basename}.fastq")
	`usearch -search_oligodb #{rev_basename}.fastq -db primer_half_fow.fasta -strand both -userout #{rev_basename}_map2.txt -userfields query+target+qstrand+diffs+tlo+thi+qlo+qhi`
end

# Read the map files and retrieve out the good sequences from it 
fow_map = File.open("#{fow_basename}_map2.txt", "r")
rev_map = File.open("#{rev_basename}_map2.txt", "r")

fow_map_hash = {}
rev_map_hash = {}
fow_retrieved_hash = {}
rev_retrieved_hash = {}

# For ones which map to the entire forward primer......
fow_map.each do |line|
	line_array = line.split("\t")
	if fow_map_hash.has_key?(line_array[0])
    		fow_map_hash[line_array[0]].push([line_array[1], line_array[2], line_array[6].chomp.to_i, line_array[7].chomp.to_i])	
  	else
    		fow_map_hash[line_array[0]] = [[line_array[1], line_array[2], line_array[6].chomp.to_i, line_array[7].chomp.to_i]]
  	end
end

fow_map_hash.each do |key, value|
	seq_length = fastq_hash[key][0].length
	#puts seq_length
	if value.length >= 2
		if value[0][0].eql?("forward") && value[0][1].eql?("+") && value[0][2].to_i < 100 && value[1][0].eql?("reverse") && value[1][1].eql?("+") && value[1][2].to_i > seq_length-100
			#puts "ONE #{value[0][3]} to #{value[1][2]}"
			fow_retrieved_hash[key] = [fastq_hash[key][0][value[0][3]..value[1][2]], fastq_hash[key][1][value[0][3]..value[1][2]]]
		elsif value[0][0].eql?("reverse") && value[0][1].eql?("+") && value[0][2].to_i > seq_length-100 && value[1][0].eql?("forward") && value[1][1].eql?("+") && value[1][2].to_i < 100
			#puts "TWO #{value[1][3]} to #{value[0][2]}"
			fow_retrieved_hash[key] = [fastq_hash[key][0][value[1][3]..value[0][2]], fastq_hash[key][1][value[1][3]..value[0][2]]]
		end
	end
end	


# For ones which map to the entire reverse primer......
rev_map.each do |line|
	line_array = line.split("\t")
	if rev_map_hash.has_key?(line_array[0])
    		rev_map_hash[line_array[0]].push([line_array[1], line_array[2], line_array[6].chomp.to_i, line_array[7].chomp.to_i])	
  	else
    		rev_map_hash[line_array[0]] = [[line_array[1], line_array[2], line_array[6].chomp.to_i, line_array[7].chomp.to_i]]
  	end
end

rev_map_hash.each do |key, value|
	seq_length = fastq_hash[key][0].length
	#puts seq_length
	if value.length >= 2
		if value[0][0].eql?("reverse") && value[0][1].eql?("+") && value[0][2].to_i > seq_length-100 && value[1][0].eql?("forward") && value[1][1].eql?("+") && value[1][2].to_i < 100
			#puts "THREE #{value[1][3]} to #{value[0][2]}"	
			rev_retrieved_hash[key] = [fastq_hash[key][0][value[1][3]..value[0][2]], fastq_hash[key][1][value[1][3]..value[0][2]]]
		elsif value[0][0].eql?("forward") && value[0][1].eql?("+") && value[0][2].to_i < 100 && value[1][0].eql?("reverse") && value[1][1].eql?("+") && value[1][2].to_i > seq_length-100
			#puts "FOUR #{value[0][3]} to #{value[1][2]}"
			rev_retrieved_hash[key] = [fastq_hash[key][0][value[0][3]..value[1][2]], fastq_hash[key][1][value[0][3]..value[1][2]]]
		end
	end
end	

puts "forward retrieved ",fow_retrieved_hash.length
puts "reverse retrieved ",rev_retrieved_hash.length

fow_retrieved_hash.each do |key, value|
	out_file_fow_retreived.puts("@#{key}")
	out_file_fow_retreived.puts(value[0])
	out_file_fow_retreived.puts("+")
	out_file_fow_retreived.puts(value[1])
end		

rev_retrieved_hash.each do |key, value|
	out_file_rev_retreived.puts("@#{key}")
	out_file_rev_retreived.puts(value[0])
	out_file_rev_retreived.puts("+")
	out_file_rev_retreived.puts(value[1])
end

