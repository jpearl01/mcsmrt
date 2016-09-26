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

out_file.puts(header_split[0..-5].join("\t")+"\tdomain\tdomain_conf\tphylum\tphylum_conf\tclass\tclass_conf\torder\torder_conf\tfamily\tfamily_conf\tgenus\tgenus_conf\tspecies\tspecies_conf\tstrain\tribo_seq_form"+"\t"+header_split[-3]+"\t"+header_split[-2]+"\t"+header_split[-1])

def parse_each_line (otublast_file, out_file)
	count = 0
	otublast_file.each do |line|
  		count+= 1
  		if count == 1
    		next
  		end
		line_split = line.split("\t")
		taxa = line_split[-4]
		capture_array = /d:([^(]+)\(([^)]+)\),p:([^(]+)\(([^)]+)\),c:([^(]+)\(([^)]+)\),o:([^(]+)\(([^)]+)\),f:([^(]+)\(([^)]+)\),g:([^(]+)\(([^)]+)\),s:([^(]+)\(([^)]+)\)/.match(line)
		#puts capture_array[1]+"\t"+capture_array[2]+"\t"+capture_array[3]+"\t"+capture_array[4]+"\t"+capture_array[5]+"\t"+capture_array[6]+"\t"+capture_array[7]+"\t"+capture_array[8]+"\t"+capture_array[9]+"\t"+capture_array[10]+"\t"+capture_array[11]+"\t"+capture_array[12]+"\t"+capture_array[13]
		#puts line_split[-3]+"\t"+capture_array[15]
		strain = ""
		ribo_type = ""
		if line_split[-3].include?("strain") and line_split[-3].include?("16S")
			#puts line_split[-3]
			strain = /(?<=strain_)(.*)(?=_16S)/.match(line_split[-3])
			#puts strain
			if line_split[-3].include?("complete")
				ribo_type = "complete"
			elsif line_split[-3].include?("partial")
				ribo_type = "partial"
			else
				ribo_type = "nil"
			end
		elsif line_split[-3].include?("strain")
			#puts line_split[-3]
			strain = /(?<=strain_)(.*)(?=_)/.match(line_split[-3])
			#puts strain
			ribo_type = "nil"
		elsif line_split[-3].include?("16S")
			#puts line_split[-3]
			strain = /(?<=_)(.*)(?=_16S)/.match(line_split[-3])
			#puts strain
			if strain.nil?
				strain = "nil"
			end
			if line_split[-3].include?("complete")
				ribo_type = "complete"
			elsif line_split[-3].include?("partial")
				ribo_type = "partial"
			else
				ribo_type = "nil"
			end
		else
			#puts line_split[-3]
			strain = /(?<=_)(.*)(?=_)/.match(line_split[-3])
			#puts strain
			ribo_type = "nil"
		end
		out_file.puts(line_split[0..-5].join("\t")+"\t"+capture_array[1]+"\t"+capture_array[2]+"\t"+capture_array[3]+"\t"+capture_array[4]+"\t"+capture_array[5]+"\t"+capture_array[6]+"\t"+capture_array[7]+"\t"+capture_array[8]+"\t"+capture_array[9]+"\t"+capture_array[10]+"\t"+capture_array[11]+"\t"+capture_array[12]+"\t"+capture_array[13]+"\t"+capture_array[14]+"\t"+strain.to_s+"\t"+ribo_type+"\t"+line_split[-3..-1].join("\t"))
	end	  
end

parse_each_line(otublast_file, out_file)
