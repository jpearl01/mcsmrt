#!/usr/bin/env ruby

require 'optimist'
require 'fileutils'

opts = Optimist::options do
        opt :fastq_folder, "Folder with fastq files to updat", type: :string, short: "-f", required: true
        opt :passes_file, "File with read header and ccs pass number: `read_header      pass_#`", type: :string, short: "-p", required: true
        opt :outfolder, "Output folder", type: :string, short: "-o", required: true
        opt :barcode_to_sample, "Rename fastq files from the barcode id to sample id", type: :string, short: "-b", required: false
end

fq_folder = opts[:fastq_folder]
passes    = opts[:passes_file]
outfolder = opts[:outfolder]
bc2samp   = opts[:barcode_to_sample] if !opts[:barcode_to_sample].nil?

FileUtils.mkdir_p outfolder unless File.exists?(outfolder)

passes_dict = {}

File.foreach(passes) do |line|
        s = line.chomp.split(/[,\t]/)
        passes_dict[s[0]]=s[1]

end
puts passes_dict.first

Dir.glob(File.join(fq_folder,'/*.f*')) do |f|
        puts f
        h = File.basename(f, ".*")
        outfile = File.open(File.join(outfolder, h+".fq"), 'w')
        File.foreach(f) do |line|
                if /^@m[0-9]+/.match(line)
                        rec = line.chomp.gsub(/^@/, '')
                        if passes_dict.has_key?(rec)
                                pass = passes_dict[rec]
                        else
                                abort("#{rec} rec did not have a number of passes")
                        end

                        line = "@#{rec};barcodelabel=#{h};ccs=#{pass};"

                end
                outfile.puts line
        end
end
