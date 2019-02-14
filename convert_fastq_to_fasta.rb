#!/usr/bin/env ruby
require 'bio'
require 'optimist'

opts = Optimist::options do
  opt :fastqfile, "fastq file with all the oriented sequences", :type => :string, :short => "-q"
  opt :fastafile, "Output fasta file", :type => :string, :short => "-a"
end

# Input - Fastq file 
# Output - Fasta file

file = Bio::FlatFile.auto(opts[:fastqfile], "r")
out_file = File.open(opts[:fastafile], "w")

file.each do |entry|
	out_file.puts(">"+entry.definition)
	out_file.puts(entry.sequence_string)
end
