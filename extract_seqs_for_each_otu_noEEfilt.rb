#!/usr/bin/env ruby
require 'bio'
require 'optimist'

##### Get the path to the directory in which the scripts exist 
script_directory = File.dirname(__FILE__)

opts = Optimist::options do
  opt :otufile, "Fasta file with all the OTUs", :type => :string, :short => "-u"
  opt :uglobfile, "File from usearch global run with all reads without EE filter", :type => :string, :short => "-p"
  opt :fqfile, "Fastq file with all the sequences", :type => :string, :short => "-a"
end

otu_file = Bio::FlatFile.auto(opts[:otufile], "r")
uglob_file = opts[:uglobfile]
fq_file = opts[:fqfile]

otu_list = []
otu_file.each do |entry|
  otu_list.push(entry.definition)
end
#puts otu_list

(0..otu_list.size-1).each do |otu|
  puts "#{otu_list[otu]}"
  otu_name = otu_list[otu].split(";")[0]
  `grep -P "#{otu_list[otu]}" #{uglob_file} | cut -f1 > split_otus/noEE_#{otu_name}_headers.txt`
  `usearch -fastx_getseqs #{fq_file} -labels split_otus/noEE_#{otu_name}_headers.txt -fastaout split_otus/noEE_#{otu_name}_eefilt_subset.fasta`
end