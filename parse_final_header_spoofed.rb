require 'bio'
require 'trollop'

opts = Trollop::options do
  opt :otublastfile, "File with OTU and blast information together", :type => :string, :short => "-h"
  opt :outfile, "Output file with headers and dereplication flag", :type => :string, :short => "-o"
end

otublast_file = File.open(opts[:otublastfile], "r")
out_file = File.open(opts[:outfile], "w")

header = File.open(opts[:otublastfile], &:readline)
header_split= header.split("\t")

out_file.puts(header_split[0..-5].join("\t")+"\tphylum\tphylum_conf\tclass\tclass_conf\torder\torder_conf\tfamily\tfamily_conf\tgenus\tgenus_conf\tspecies\tspecies_conf\tstrain\tstrain_conf"+"\t"+header_split[-3]+"\t"+header_split[-2]+"\t"+header_split[-1])

def parse_each_line (otublast_file, out_file)
	count = 0
	otublast_file.each do |line|
  		count+= 1
  		if count == 1
    		next
  		end
		line_split = line.split("\t")
		taxa = line_split[-4]
		#puts taxa
		capture_array = /d:([^(]+)\(([^)]+)\),p:([^(]+)\(([^)]+)\),c:([^(]+)\(([^)]+)\),o:([^(]+)\(([^)]+)\),f:([^(]+)\(([^)]+)\),g:([^(]+)\(([^)]+)\),s:([^_]+_[^_]+)(.*)\(([^)]+)/.match(line)
		#puts capture_array[1]+"\t"+capture_array[2]+"\t"+capture_array[3]+"\t"+capture_array[4]+"\t"+capture_array[5]+"\t"+capture_array[6]+"\t"+capture_array[7]+"\t"+capture_array[8]+"\t"+capture_array[9]+"\t"+capture_array[10]+"\t"+capture_array[11]+"\t"+capture_array[12]+"\t"+capture_array[13]
		#puts capture_array[14]+"\t"+capture_array[15]
		out_file.puts(line_split[0..-5].join("\t")+"\t"+capture_array[1]+"\t"+capture_array[2]+"\t"+capture_array[3]+"\t"+capture_array[4]+"\t"+capture_array[5]+"\t"+capture_array[6]+"\t"+capture_array[7]+"\t"+capture_array[8]+"\t"+capture_array[9]+"\t"+capture_array[10]+"\t"+capture_array[11]+"\t"+capture_array[12]+"\t"+capture_array[14]+"\t"+capture_array[15]+"\t"+line_split[-3..-1].join("\t"))	
	end	  
end

parse_each_line(otublast_file, out_file)
