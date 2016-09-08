require 'bio'

sample = ARGV[0]

class Filtered_steps_2
  attr_accessor :ee_count, :derep_count, :chimera_in_clustering, :chimera_in_uchime, :ambiguous_in_uchime, :num_of_otus
end

# File mapping between the OTU table and the utax table                                                                                                     
table_utax_map_file = File.open("#{sample}_table_utax_map.txt", 'w')
report_file = File.open("#{sample}_report.txt", 'w')

#####Initialize log variable                                                                                                                                
filt_2 = Filtered_steps_2.new
filt_2.ee_count=0
filt_2.derep_count=0
filt_2.chimera_in_clustering=0
filt_2.chimera_in_uchime=0
filt_2.ambiguous_in_uchime=0

##### Dealing with the file after expected error filtering  
oriented_filtered_open = Bio::FlatFile.auto("#{sample}_oriented_filtered.fastq")
count_ee = oriented_filtered_open.to_a.count
filt_2.ee_count = count_ee

##### Dealing with the dereplicated file 
derep_open =  Bio::FlatFile.auto("#{sample}_oriented_filtered_derep.fasta")                                                                                                     
count_derep = derep_open.to_a.count
filt_2.derep_count = count_derep

##### Dealing with the up file from clustering step                                                                                                          
count_up = `grep "chimera" #{sample}_uparse_out.up -c`
filt_2.chimera_in_clustering = count_up.to_i

##### Dealing with the chimeras file
if File.zero?("#{sample}_OTU_chimeras.fasta")                                                                                                               
  	filt_2.chimera_in_uchime = 0
else
	chimera_open = Bio::FlatFile.auto("#{sample}_OTU_chimeras.fasta")
	count_chimera = chimera_open.to_a.count
  	filt_2.chimera_in_uchime = count_chimera
end 

##### Dealing with uchime file                                                                                                                               
count_uchime = `grep "?" #{sample}_OTU_uchime.txt -c`
uchime_file = File.open("#{sample}_OTU_uchime.txt", 'r')
otu_size = []
query = []
parent = []
id = []
uchime_file.each do |line|
  if line.include?("?")
    uchime_line = line.split("\t")
    otu_size.push(uchime_line[1].split(";")[1])
    query.push(uchime_line[2].split("|")[1])
    parent.push(uchime_line[3].split("|")[1])
    id.push(uchime_line[4])
  end  
end
filt_2.ambiguous_in_uchime = count_uchime.to_i

##### Dealing with OTU table and utax file.... Mapping the two!
otus = `wc -l #{sample}_Table.txt`
filt_2.num_of_otus = otus.split(" ")[0]                                                                                              
otu_table_file = File.open("#{sample}_Table.txt", "r")
otu_table_hash = {}
otu_table_file.each do |line|
  line_array = line.split("\t")
  otu_table_hash[line_array[0].chomp] = line_array[1..-1]
end
#puts otu_table_hash.sort                                                                                                                                         

utax_table_file = File.open("#{sample}_reads.utax", "r")
utax_table_hash = {}
utax_table_hash["#OTU ID"] = "Assignment"
size_hash = {}
size_hash["#OTU ID"] = "size"
utax_table_file.each do |line|
  line_array = line.split("\t")
  utax_table_hash[line_array[0].split(";")[0].chomp] = line_array[1].chomp
  size_hash[line_array[0].split(";")[0].chomp] = line_array[0].split(";")[1].split("=")[1].chomp
end
#puts utax_table_hash.sort
#puts size_hash                                                                                                                                      

otu_table_hash.sort.each do |key, value|
  	value_to_print = value.join("\t").chomp
  	table_utax_map_file.puts("#{key}\t#{size_hash[key]}\t#{value_to_print}\t#{utax_table_hash[key]}")
end

##### Printing the final report                                                                                                                              
report_file.puts("Sample\tEE_count\tDereplication_count\tChimeras_from_clustering\tChimeras_from_uchime\tAmbiguous_otus_in_uchime\tOTU_size\tNum_of_OTUs")

if otu_size.size == 0
  report_file.puts("#{sample}\t#{filt_2.ee_count}\t#{filt_2.derep_count}\t#{filt_2.chimera_in_clustering}\t#{filt_2.chimera_in_uchime}\t#{filt_2.ambiguous_in_uchime}\tNA\t#{filt_2.num_of_otus}")
else
  report_file.puts("#{sample}\t#{filt_2.ee_count}\t#{filt_2.derep_count}\t#{filt_2.chimera_in_clustering}\t#{filt_2.chimera_in_uchime}\t#{filt_2.ambiguous_in_uchime}\t#{otu_size*","}\t#{filt_2.num_of_otus}")
end


