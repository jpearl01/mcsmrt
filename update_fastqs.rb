#!/usr/bin/env ruby

require 'optimist'


opts = Optimist::options do
 opt :foldername, "Folder with fastq files, named like 'lima_output.lbc30--lbc30.bam.Q20.fastq'", type: :string, short: "-f"
end

folder_path = opts[:foldername]

Dir.glob(File.join(opts[:foldername],"*")) do |f|
	puts "fixing #{f}"
	new_name = /[^\.]+\.+([^-]+)[^\.]+\.bam.+fastq/.match(f)[1]
	File.rename(f, File.join(folder_path,  new_name + ".fq"))
end