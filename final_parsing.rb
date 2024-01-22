require 'optimist'

opts = Optimist::options do
  opt :usearchfile, "File with usearch local information.", :type => :string, :short => "-b"
  opt :otuutaxfile, "File with OTU table and utax information.", :type => :string, :short => "-u", required: true
  opt :ncbiclusteredfile, "File with info about DB clustered", :type => :string, :short => "-n"
  opt :outfile, "Output file with OTU names, OTU counts, lineage and blast results.", :type => :string, :short => "-o", required: true
end

usearch_file = File.open(opts[:usearchfile], "r") unless opts[:usearchfile].nil?
otu_utax_file = File.open(opts[:otuutaxfile], "r")
ncbi_clustered_file = File.open(opts[:ncbiclusteredfile], "r") unless opts[:ncbiclusteredfile].nil?
out_file = File.open(opts[:outfile], "w")

usearch_file_hash = {}
if !usearch_file.nil?
	usearch_file.each do |line|
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
		usearch_file_hash[key] = [query, strain, complete, percent_identity, alignment_length]
	end
	#puts usearch_file_hash
end

ncbi_clustered_hash = {}
if !ncbi_clustered_file.nil?
	ncbi_clustered_file.each_with_index do |line, index|
		if index == 0
			next
		else
			line_split = line.split("\t")
			species = line_split[12]
			num_of_species = line_split[-3]
			avg_otu_id = line_split[1]
			num_of_avg_otu_with_same_species = line_split[-1]
			#puts species
			if ncbi_clustered_hash.has_key?(species)
				ncbi_clustered_hash[species][2] += 1
			else
				ncbi_clustered_hash[species] = [avg_otu_id, num_of_species, 1, num_of_avg_otu_with_same_species]
			end
		end 
	end
end

#puts ncbi_clustered_hash["Microbacter_margulisiae"]
#puts ncbi_clustered_hash

header = File.open(otu_utax_file, &:readline)
header_split = header.split("\t")[0..-2].join("\t")
out_file.puts(header_split+"\tdomain\tdomain_conf\tphylum\tphylum_conf\tclass\tclass_conf\torder\torder_conf\tfamily\tfamily_conf\tgenus\tgenus_conf\tspecies\tspecies_conf\tusearch_query\tusearch_strain\tusearch_16s_completeness\tusearch_percent_identity\tusearch_alignment_length\tncbi_avglinkage_otu_id\tnum_of_sp_in_cluster\tnum_of_sp_in_db\tnum_of_otus_with_same_sp")

otu_utax_file.each_with_index do |line, index|
	if index == 0
		#puts line
		next
	else
	  line_split = line.split("\t")
    key = line_split[0]
    capture_array = /d:([^(]+)\(([^)]+)\),p:([^(]+)\(([^)]+)\),c:([^(]+)\(([^)]+)\),o:([^(]+)\(([^)]+)\),f:([^(]+)\(([^)]+)\),g:([^(]+)\(([^)]+)\),s:([^(]+)\(([^)]+)\)/.match(line)
    #puts capture_array
    puts capture_array[1]+"\t"+capture_array[2]+"\t"+capture_array[3]+"\t"+capture_array[4]+"\t"+capture_array[5]+"\t"+capture_array[6]+"\t"+capture_array[7]+"\t"+capture_array[8]+"\t"+capture_array[9]+"\t"+capture_array[10]+"\t"+capture_array[11]+"\t"+capture_array[12]+"\t"+capture_array[13]
    #next unless usearch_file_hash.has_key?(key)
    if ncbi_clustered_hash.has_key?(capture_array[13]) and usearch_file_hash.has_key?(key)
     	out_file.puts(line_split[0..-2].join("\t")+"\t"+capture_array[1..-1].join("\t")+"\t"+usearch_file_hash[key][0..-1].join("\t")+"\t"+ncbi_clustered_hash[capture_array[13]].join("\t"))
		elsif usearch_file_hash.has_key?(key)
			out_file.puts(line_split[0..-2].join("\t")+"\t"+capture_array[1..-1].join("\t")+"\t"+usearch_file_hash[key][0..-1].join("\t")+"\tNA\tNA\tNA\tNA")
		elsif ncbi_clustered_hash.has_key?(capture_array[13]) 
			out_file.puts(line_split[0..-2].join("\t")+"\t"+capture_array[1..-1].join("\t")+"\tNA\tNA\tNA\tNA\tNA\t"+ncbi_clustered_hash[capture_array[13]].join("\t"))
		else
			out_file.puts(line_split[0..-2].join("\t")+"\t"+capture_array[1..-1].join("\t")+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")
		end
	end	
end
