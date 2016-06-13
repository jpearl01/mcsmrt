require 'trollop'

opts = Trollop::options do
	opt :samplecountsfile, "File created from the main pipeline with the counts at each step", :type => :string, :short => "-s"
	opt :outfile, "Filtered counts with the leveled output at each step (reducing)", :type => :string, :short => "-o"
end
#puts opts

sample_counts_file = File.open(opts[:samplecountsfile], "r")
out_file = File.open(opts[:outfile], "w")

out_file.puts("Sample\tOriginal_count\tSize_filt\tHuman_map\tUSERACH_primer_match\tCustom_primer_match_Orienting")

count = 0
sample_counts_file.each do |line|
	count+= 1
	if count == 1
		next
	else
		line_split = line.split("\t")
		sample = line_split[0]
		og_count = line_split[1]
		size_filt = line_split[1].to_i-(line_split[2].to_i+line_split[3].to_i)
		human_map = size_filt-line_split[4].to_i	
		not_primer_mapped_usearch = human_map-line_split[8].to_i
		lost_in_custom_primer_map = (line_split[5].to_i+line_split[6].to_i)-line_split[10].to_i
		custom_primer_map = not_primer_mapped_usearch-lost_in_custom_primer_map
		out_file.puts(sample+"\t"+og_count+"\t"+size_filt.to_s+"\t"+human_map.to_s+"\t"+not_primer_mapped_usearch.to_s+"\t"+custom_primer_map.to_s)
	end
end
	
out_file.close
