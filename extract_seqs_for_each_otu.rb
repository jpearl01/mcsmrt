#!/usr/bin/env ruby
require 'bio'
require 'optimist'

opts = Optimist::options do
  opt :otufile, "Fasta file with all the OTUs", :type => :string, :short => "-u"
  opt :upfile, "File from uparse", :type => :string, :short => "-p"
  opt :fafile, "Fasta file with all the dereplicated sequences", :type => :string, :short => "-a"
end

otu_file = Bio::FlatFile.auto(opts[:otufile], "r")
up_file = opts[:upfile]
fa_file = opts[:fafile]

otu_list = []
otu_file.each do |entry|
  otu_list.push(entry.definition.split(";")[0])
end
#puts otu_list

(0..otu_list.size-1).each do |otu|
	otu_name = otu_list[otu].split(";")[0]
  	puts "#{otu_list[otu]}"
  	`grep -P "#{otu_list[otu]}$" #{up_file} | cut -f1 > split_otus/#{otu_name}_headers.tsv`
  	`usearch -fastx_getseqs #{fa_file} -labels split_otus/#{otu_name}_headers.tsv -fastaout split_otus/#{otu_name}_eefilt_subset.fasta`
end