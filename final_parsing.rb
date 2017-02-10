require 'trollop'

opts = Trollop::options do
  opt :blastfile, "File with blast information.", :type => :string, :short => "-b"
  opt :otuutaxfile, "File with OTU table and utax information.", :type => :string, :short => "-u"
  opt :outfile, "Output file with OTU names, OTU counts, lineage and blast results.", :type => :string, :short => "-o"
end

blast_file = File.open(opts[:blastfile], "r")
otu_utax_file = File.open(opts[:otuutaxfile], "r")
out_file = File.open(opts[:outfile], "w")

blast_file_hash = {}
blast_file.each do |line|
	line_split = line.split("\t")
	#puts line_split[0]
	# get the otu name which will be the key of the hash
	key = line_split[0].split(";")[0]
	# Check if strain info is there in query, also get the query name! 
	if line.include?("strain=")
		strain = line_split[1].split(";")[2].split("=")[1]
		query = line_split[1].split(";")[1].split("=")[1]
	else
		strain = "NA"
		query = line_split[1]
	end
	# Check if 16s completeness is there in query
	if line.include?("complete=")
		complete = line_split[1].split(";")[3].split("=")[1]
	else
		complete = "NA"
	end
	# Get the percent identity and alignment length
	percent_identity = line_split[2]
	alignment_length = line_split[3]

	# populate the hash 
	blast_file_hash[key] = [query, strain, complete, percent_identity, alignment_length]
end
#puts blast_file_hash

header = File.open(otu_utax_file, &:readline)
header_split = header.split("\t")[0..-2].join("\t")
out_file.puts(header_split+"\tdomain\tdomain_conf\tphylum\tphylum_conf\tclass\tclass_conf\torder\torder_conf\tfamily\tfamily_conf\tgenus\tgenus_conf\tspecies\tspecies_conf\tblast_query\tblast_strain\tblast_16s_completeness\tblast_percent_identity\tblast_alignment_length")

otu_utax_file.each_with_index do |line, index|
	if index == 0
		next
	else
		line_split = line.split("\t")
		key = line_split[0]
		capture_array = /d:([^(]+)\(([^)]+)\),p:([^(]+)\(([^)]+)\),c:([^(]+)\(([^)]+)\),o:([^(]+)\(([^)]+)\),f:([^(]+)\(([^)]+)\),g:([^(]+)\(([^)]+)\),s:([^(]+)\(([^)]+)\)/.match(line)
		#puts capture_array[1]+"\t"+capture_array[2]+"\t"+capture_array[3]+"\t"+capture_array[4]+"\t"+capture_array[5]+"\t"+capture_array[6]+"\t"+capture_array[7]+"\t"+capture_array[8]+"\t"+capture_array[9]+"\t"+capture_array[10]+"\t"+capture_array[11]+"\t"+capture_array[12]+"\t"+capture_array[13]
		out_file.puts(line_split[0..-2].join("\t")+"\t"+capture_array[1..-1].join("\t")+"\t"+blast_file_hash[key][0..-1].join("\t"))
	end	
end