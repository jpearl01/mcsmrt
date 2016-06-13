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

out_file.puts(header_split[0..-5].join("\t")+"\tdomain\tdomain_conf\tphylum\tphylum_conf\tclass\tclass_conf\torder\torder_conf\tfamily\tfamily_conf\tgenus\tgenus_conf\tspecies\tstrain\tseq_form"+"\t"+header_split[-3]+"\t"+header_split[-2]+"\t"+header_split[-1])

count = 0
otublast_file.each do |line|
  count+= 1
  if count == 1
    next
  end  
  line_split = line.split("\t")
  #puts line_split[-4]
  # Dealing with the assignments and taxonomy 
  if line_split[-4].include?("d:")
    domain = /d\:(.+?),/.match(line).to_s.split(":")[1].split("(")[0]
    domain_conf = /d\:(.+?),/.match(line).to_s.split(":")[1].split("(")[1].tr("),", "")
  else
    domain = "nil"
    domain_conf = "nil"
  end
  if line_split[-4].include?("p:")
    phylum = /p\:(.+?),/.match(line).to_s.split(":")[1].split("(")[0]
    phylum_conf = /p\:(.+?),/.match(line).to_s.split(":")[1].split("(")[1].tr("),", "") 
  else
    phylum = "nil"
    phylum_conf = "nil"
  end
  if line_split[-4].include?("c:")
    clas = /c\:(.+?),/.match(line).to_s.split(":")[1].split("(")[0]
    clas_conf = /c\:(.+?),/.match(line).to_s.split(":")[1].split("(")[1].tr("),", "")
  else
    clas = "nil"
    clas_conf = "nil"
  end
  if line_split[-4].include?("o:")
    order = /o\:(.+?),/.match(line).to_s.split(":")[1].split("(")[0]
    order_conf = /o\:(.+?),/.match(line).to_s.split(":")[1].split("(")[1].tr("),", "")
  else
    order = "nil"
    order_conf = "nil"
  end
  if line_split[-4].include?("f:")
    family = /f\:(.+?),/.match(line).to_s.split(":")[1].split("(")[0]
    family_conf = /f\:(.+?),/.match(line).to_s.split(":")[1].split("(")[1].tr("),", "")
  else
    family = "nil"
    family_conf = "nil"
  end
  if line_split[-4].include?("g:")
    genus = /g\:(.+?),/.match(line).to_s.split(":")[1].split("(")[0]
    genus_conf = /g\:(.+?),/.match(line).to_s.split(":")[1].split("(")[1].tr("),", "")
  else
    genus = "nil"
    genus_conf = "nil"
  end
  if line_split[-4].include?("s:")
    species = /s\:(.+?);/.match(line).to_s.split(":")[1].split("(")[0]
    species_to_print = /s\:(.+?);/.match(line).to_s.split(":")[1].split("(")[0].split("_")[0]+"_"+/s\:(.+?);/.match(line).to_s.split(":")[1].split("(")[0].split("_")[1]
    if species.include?("strain")
      strain = /strain_(.*)_/.match(species).to_s.split("_")[1]+"_"+/strain_(.*)_/.match(species).to_s.split("_")[2]
    else
      strain = "nil"
    end  
    if species.include?("complete")
      seq_form = "complete"
    elsif species.include?("partial")
      seq_form = "partial"
    else
      seq_form = "nil"
    end  
  else
    species = "nil"
  end
  out_file.puts(line_split[0..-5].join("\t")+"\t"+domain+"\t"+domain_conf+"\t"+phylum+"\t"+phylum_conf+"\t"+clas+"\t"+clas_conf+"\t"+order+"\t"+order_conf+"\t"+family+"\t"+family_conf+"\t"+genus+"\t"+genus_conf+"\t"+species_to_print+"\t"+strain+"\t"+seq_form+"\t"+line_split[-3]+"\t"+line_split[-2]+"\t"+line_split[-1])
end  
