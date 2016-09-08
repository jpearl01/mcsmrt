require 'bio'
require 'trollop'

opts = Trollop::options do
  opt :otu_file, "File created by mapping between the OTU table and utax", :type => :string, :short => "-u"
  opt :blast_file, "Blast results for each otu", :type => :string, :short => "-b"
  opt :outfile, "Output file to write the OTU sequences with the ones belonging to the ?", :type => :string, :short => "-o"
end

otu_file = File.open(opts[:otu_file], "r")
blast_file = File.open(opts[:blast_file], "r")
out_file = File.open(opts[:outfile], "w")

# Created a hash of the otu table with the otu name as the key and the entire lines as the value                                                            
otu_hash = {}
count = 0
header = ""
otu_file.each do |line|
  if count == 0
    header = line.chomp
    count += 1
    next
  else
    line_split = line.split("\t")
    count += 1
    otu_key = line_split[0]
    otu_hash[otu_key] = line.chomp
  end
end
#puts otu_hash      

# Creating a hash of the blast table with otu name as the key and identity and alignment length as the value
blast_hash = {}
blast_file.each do |line|
  line_split = line.split("\t")
  otu_key = line_split[0].split(";")[0]
  blast_hash[otu_key] = "#{line_split[1]}\t#{line_split[2]}\t#{line_split[3]}"
end
#puts blast_hash

# Creating a final hash with mapped otus and blast identities
new_hash = {}
otu_hash.each do |key, value|
  key_split = key.split(",")
  if key_split.length > 1
    otu_names = []
    size = []
    (0..key_split.length-1).each do |each_otu|
      #puts key_split[each_otu]
      otu_names.push(key_split[each_otu].split(";")[0])
      size.push((key_split[each_otu].split(";")[1].split("=")[1]).to_i) 
    end  
    #print otu_names, size
    max_size_index =  size.index(size.max)
    max_size_otu = otu_names[max_size_index]
    blast_key = "#{max_size_otu};size=#{size.max};"
    #puts blast_key
    #puts blast_hash[key_split[each_otu]]
    #puts blast_hash.keys
    if blast_hash.has_key?(blast_key)
      #print key, blast_hash[blast_key]
      new_hash[key] = "#{otu_hash[key]}\t#{blast_hash[blast_key]}"
    end
  else
    #puts blast_hash[key]
    new_hash[key] = "#{otu_hash[key]}\t#{blast_hash[key]}" 
  end  
end
#puts new_hash

out_file.puts("#{header.chomp}Assignment\tQuery\tIdentity\tAlignment_Length")
new_hash.each do |key, value|
  out_file.puts(value)
end
