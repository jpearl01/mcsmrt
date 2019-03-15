#!/usr/bin/env ruby

all_bc_info = File.open("pre_all_reads_info.tsv")
all_bc_utax = File.open("all_bc_reads.utax")
outfile = File.open("pre_all_bc_reads_info_with_taxa.tsv", "w")

utax_header = "Domain\tDomain_conf\tPhylum\tPhylum_conf\tClass\tClass_conf\tOrder\tOrder_conf\tFamily\tFamily_conf\tGenus\tGenus_conf\tSpecies\tSpecies_conf"
utax_hash = {}
all_bc_utax.each do |line|
  line_split = line.split("\t")
  capture_array = /d:([^(]+)\(([^)]+)\),p:([^(]+)\(([^)]+)\),c:([^(]+)\(([^)]+)\),o:([^(]+)\(([^)]+)\),f:([^(]+)\(([^)]+)\),g:([^(]+)\(([^)]+)\),s:([^(]+)\(([^)]+)\)/.match(line_split[1])
  key_utax = line_split[0].split(";")[0]
  if capture_array.nil?
    next
  else
    utax_hash[key_utax] = capture_array[1..-1].join("\t")
  end
end

all_reads_hash = {}
info_header = ""
all_bc_info.each_with_index do |line, index|
  if index == 0
    info_header = line.chomp
  else
    line_split = line.split("\t")
    all_reads_hash[line_split[0]] = line_split[1..-1].join("\t").chomp
  end
end

outfile.puts("#{info_header}\t#{utax_header}")
all_reads_hash.each do |key, value|
  if utax_hash.key?(key)
    outfile.puts("#{key}\t#{all_reads_hash[key]}\t#{utax_hash[key]}")
  else
    outfile.puts("#{key}\t#{all_reads_hash[key]}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")
  end
end