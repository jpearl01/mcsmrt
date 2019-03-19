#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
require 'bio'
require 'optimist'
require 'colorize'

# USAGE: ruby /data/shared/homes/archana/projects/mcsmrt/mcsmrt.rb -a -f reads/ -d 32 -e 1 -s 5 -x 2000 -n 500 -c /data/shared/homes/archana/projects/rdp_gold.fa -t /data/shared/homes/archana/projects/lineanator/16sMicrobial_ncbi_lineage_reference_database.udb -l /data/shared/homes/archana/projects/lineanator/16sMicrobial_ncbi_lineage.fasta -g /data/shared/homes/archana/projects/human_g1k_v37.fasta -p /data/shared/homes/archana/projects/primers.fasta -b /data/shared/homes/archana/projects//mcsmrt/ncbi_clustered_table.tsv -v


opts = Optimist::options do
  opt :foldername, "Folder with the demultiplexed files for clustering", :type => :string, :short => "-f"
  opt :samplelist, "File with a subset of fastq files to be used (tab separated sample_name/filepath)", :type => :string, :short => "-i"
  opt :threads, "Number of threads you can allot for running this process", :type => :int, :short => "-d", :default => 1
  opt :eevalue, "Expected error value greater than which reads will be filtered out", :type => :float, :short => "-e", :default => 1.0
  opt :trimming, "Do you want to trim your sequences? Answer in yes or no", :type => :string, :short => "-m", :default => "yes"
  opt :ccsvalue, "CCS passes lesser than which reads will be filtered out", :type => :int, :short => "-s", :default => 5
  opt :lengthmax, "Maximum length above which reads will be filtered out", :type => :int, :short => "-x", :default => 2000
  opt :lengthmin, "Manimum length below which reads will be filtered out", :type => :int, :short => "-n", :default => 500
  opt :uchimedbfile, "Path to database file for the uchime command", :type => :string, :short => "-c", required: true
  opt :utaxdbfile, "Path to database file for the utax command", :type => :string, :short => "-t", required: true
  opt :lineagefastafile, "Path to FASTA file with lineage info for the ublast command", :type => :string, :short => "-l"
  opt :host_db, "Path to fasta file of host genome", :type => :string, :short => "-g"
  opt :primerfile, "Path to fasta file with the primer sequences", :type => :string, :short => "-p", required: true
  opt :ncbiclusteredfile, "Path to a file with database clustering information", :type => :string, :short => "-b"
  opt :verbose, "Use -v if you want all the intermediate files, else the 6 important results file will be kept and the rest will be deleted", :short => "-v"
  opt :splitotu, "Do you want to split OTUs into individual fasta files? Answer in yes or no", :type => :string, :short => "-o", :default => "no"
  opt :splitotumethod, "If you chose yes to split OTUs, do it before or after EE filtering? Answer in before or after", :type => :string, :short => "-j", :default => "after"
end 

##### Assigning variables to the input and make sure we got all the inputs
if opts[:samplelist] != nil  
  sample_list = opts[:samplelist]
  all_files = nil
else
  all_files = "*"
  sample_list = nil
end              


if opts[:verbose] == true     
  verbose = true   
else
  verbose = nil   
end  

File.exists?(opts[:foldername])          ? folder_name = opts[:foldername]              : abort("Must supply the name of the folder with demultiplexed fastq files'-f'")
File.exists?(opts[:uchimedbfile])        ? uchime_db_file = opts[:uchimedbfile]         : abort("Must supply an existing 'uchime database file' e.g. rdpgold.udb with '-c'")
File.exists?(opts[:utaxdbfile])          ? utax_db_file = opts[:utaxdbfile]             : abort("Must supply an existing 'utax database file' e.g. 16s_ncbi.udb with '-t'")
File.exists?(opts[:primerfile])          ? primer_file = opts[:primerfile]              : abort("Must supply an existing fasta of the primer sequences e.g primer_seqs.fa with '-p'")
File.exists?(opts[:ncbiclusteredfile])   ? ncbi_clust_file = opts[:ncbiclusteredfile]   : abort("Must supply an existing file with database clustering information with '-b'")

lineage_fasta_file = nil
if !opts[:lineagefastafile].nil?
  File.exists?(opts[:lineagefastafile])    ? lineage_fasta_file = opts[:lineagefastafile] : abort("Must supply an existing 'lineage fasta file' e.g. ncbi_lineage.fasta (for blast) with '-l'")
end

human_db = nil
if !opts[:host_db].nil? 
  if File.exists?(opts[:host_db])  
    human_db = opts[:host_db]      
  else 
    abort("Cannot find host genome file e.g. human_g1k.fasta with '-g'")
  end
end



thread = opts[:threads].to_i
ee = opts[:eevalue].to_f 
trim_req = opts[:trimming].to_s
ccs = opts[:ccsvalue].to_i 
length_max = opts[:lengthmax].to_i 
length_min = opts[:lengthmin].to_i
split_otus = opts[:splitotu].to_s
split_otu_method = opts[:splitotumethod].to_s

##### Get the path to the directory in which the scripts exist 
script_directory = File.dirname(__FILE__)

##### Class that stores information about each record from the reads file
class Read_sequence
  attr_accessor :read_name, :basename, :ccs, :barcode, :sample, :ee_pretrim, :ee_posttrim, :length_pretrim, :length_posttrim, :host_map, :f_primer_matches, :r_primer_matches, :f_primer_start, :f_primer_end, :r_primer_start, :r_primer_end, :read_orientation, :primer_note, :num_of_primerhits

  def initialize
    @read_name         = "NA" 
    @basename          = "NA"
    @ccs               = -1
    @barcode           = "NA" 
    @sample            = "NA" 
    @ee_pretrim        = -1
    @ee_posttrim       = -1 
    @length_pretrim    = -1
    @length_posttrim   = -1
    @host_map          = "NA"
    @f_primer_matches  = false
    @r_primer_matches  = false
    @f_primer_start    = "NA"
    @f_primer_end      = "NA"
    @r_primer_start    = "NA"
    @r_primer_end      = "NA"
    @read_orientation  = "NA"
    @primer_note       = "NA"
    @num_of_primerhits = "NA"
  end
end

##### Class that stores information about each sample and the number of reads which pass each step
class Report
  attr_accessor :total, :ee_filt, :size_filt, :host_filt, :primer_filt, :derep_filt, :chimera_filt

  def initialize
    @total = 0
    @ee_filt = 0
    @size_filt = 0
    @host_filt = 0
    @primer_filt = 0
    @derep_filt = 0
    @chimera_filt = 0
  end
end


##################### METHODS #######################

##### Method to concatenate files to create one file for clustering
def concat_files (folder_name, all_files, sample_list)
  #puts all_files
  if all_files.nil?
    sample_list_file = File.open(sample_list)
    list = []
    sample_list_file.each do |line|
      list.push("#{folder_name}/"+line.strip)
    end
    #puts list
    `cat #{list.join(" ")} > pre_demultiplexed_ccsfilt.fq`

  else
    `cat #{folder_name}/#{all_files} > pre_demultiplexed_ccsfilt.fq` 
  end

  abort("!!!!The file with all the reads required for clustering does not exist!!!!") if !File.exists?("pre_demultiplexed_ccsfilt.fq")
end

##### Method to write reads in fastq format
#fh = File handle, header = read header (string), sequence = read sequence, quality = phred quality scores
def write_to_fastq (fh, header, sequence, quality)
=begin
  fh       = file handle
  header   = string
  sequence = string
  quality  = array
=end
  fh.write('@' + header + "\n")
  fh.write(sequence)
  fh.write("\n+\n")
  fh.write(quality + "\n")
end

##### Method whcih takes an fq file as argument and returns a hash with the read name and ee
def get_ee_from_fq_file (file_basename, ee, suffix)
	`usearch -fastq_filter #{file_basename}.fq -fastqout #{suffix} -fastq_maxee 20000 -fastq_qmax 127 -fastq_qmaxout 127 -fastq_eeout -sample all`
	
  abort("!!!!Expected error filtering with usearch failed!!!!") if File.zero?("#{suffix}")

	# Hash that is returned from this method (read name - key, ee - value)
	ee_hash = {}
	
	# Open the file from fastq_filter command and process it
	ee_filt_file = Bio::FlatFile.auto("#{suffix}")
	ee_filt_file.each do |entry|
		entry_def_split = entry.definition.split(";")
		read_name = entry_def_split[0].split("@")[0]
		ee = entry.definition.match(/ee=(.*);/)[1]
		ee_hash[read_name] = ee.to_f
	end

  #abort("!!!!Dictionary from expected error filtering with usearch is empty!!!!") if ee_hash.empty?("#{file_basename}_#{suffix}")
	
	return ee_hash
end

##### Mapping reads to the human genome
def map_to_human_genome (file_basename, human_db, thread)
  #align all reads to the human genome                                                                                                                   
  `bwa mem -t #{thread} #{human_db} #{file_basename}.fq > pre_map_to_host.sam`

  # Check to make sure bwa worked, which means that the index files also exist
  abort("!!!!Index files are missing, run the bwa index command to index the reference genome and store them in the same directory as the reference genome!!!!") if File.zero?("pre_map_to_host.sam")
  
  #sambamba converts sam to bam format                                                                                                                   
  `sambamba view -S -f bam pre_map_to_host.sam -o pre_map_to_host.bam`
  
  #Sort the bam file                                                                                                                                     
  #`sambamba sort -t#{thread} -o filt_non_host_sorted.bam filt_non_host.bam`
  
  #filter the bam for only ‘not unmapped’ reads -> reads that are mapped                                                                                 
  `sambamba view -F 'not unmapped' pre_map_to_host.bam > pre_host_mapped.txt`
  mapped_count = `cut -d ';' -f1 pre_host_mapped.txt | sort | uniq | wc -l`
  mapped_string = `cut -d ';' -f1 pre_host_mapped.txt`
  
  #filter reads out for ‘unmapped’ -> we would use these for pipeline                                                                             
  `sambamba view -F 'unmapped' pre_map_to_host.bam > pre_filt_non_host.txt`
  
  #convert the sam file to fastq                                                                                                                         
  `grep -v ^@ pre_filt_non_host.txt | awk '{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}' > pre_filt_non_host.fq`
  
 	return mapped_count, mapped_string
end

##### Method for primer matching 
def primer_match (script_directory, file_basename, primer_file, thread)
	# Run the usearch command for primer matching
  puts "usearch -usearch_local #{file_basename}.fq -db #{primer_file} -id 0.8 -threads #{thread} -userout pre_primer_map.tsv -userfields query+target+qstrand+tlo+thi+qlo+qhi+pairs+gaps+mism+diffs -strand both -gapopen 1.0 -gapext 0.5 -lopen 1.0 -lext 0.5"
	`usearch -usearch_local #{file_basename}.fq -db #{primer_file} -id 0.8 -threads #{thread} -userout pre_primer_map.tsv -userfields query+target+qstrand+tlo+thi+qlo+qhi+pairs+gaps+mism+diffs -strand both -gapopen 1.0 -gapext 0.5 -lopen 1.0 -lext 0.5`

  # Check to see if sequences passed primer matching, i.e., if no read has a hit for primers, the output file from the previous step will be empty!
  abort("!!!!None of the reads mapped to the primers, check your FASTA file which has the primers!!!!") if File.zero?("pre_primer_map.tsv")

  # Run the script which parses the primer matching output
  `ruby #{script_directory}/primer_parse.rb -p pre_primer_map.tsv -o pre_primer_map_info.tsv` 
  	
  # Open the file with parsed primer matching results
  "pre_primer_map_info.tsv".nil? ==false  ? primer_matching_parsed_file = File.open("pre_primer_map_info.tsv") : abort("Primer mapping parsed file was not created from primer_parse.rb")

  return_hash = {}
  primer_matching_parsed_file.each_with_index do |line, index|
    if index == 0
      next
    else
      line_split = line.chomp.split("\t")
      key = line_split[0]
		  return_hash[key] = Array.new(line_split[1..-1].length)
		  if line_split[1] == "true"
        return_hash[key][0] = true
		  else
  		  return_hash[key][0] = false
      end
		  if line_split[2] == "true"
        return_hash[key][1] = true
      else
        return_hash[key][1] = false
      end
      return_hash[key][2..-1] = line_split[3..-1]
    end
  end
  return return_hash
end

##### Method to trim and orient the sequences
def trim_and_orient (f_primer_matches, r_primer_matches, f_primer_start, f_primer_end, r_primer_start, r_primer_end, read_orientation, seq, qual)
  # Variable which will have the seq and qual stings oriented and trimmed
  seq_trimmed = ""
  qual_trimmed = ""

  # For cases when both primers matched
  if f_primer_matches == true and r_primer_matches == true
    if read_orientation == "+"
      seq_trimmed = seq[f_primer_end..r_primer_start]
      qual_trimmed = qual[f_primer_end..r_primer_start]
    elsif read_orientation == "-"
      seq = Bio::Sequence::NA.new(seq)
      revcom_seq = seq.complement.upcase
      rev_qual = qual.reverse
      seq_trimmed = revcom_seq[r_primer_end..f_primer_start]
      qual_trimmed = rev_qual[r_primer_end..f_primer_start]
    end
  end

  if seq_trimmed != ""
    ee = 0.00
    seq_trimmed_bio = Bio::Fastq.new("@name\n#{seq_trimmed}\n+\n#{qual_trimmed}")
    quality_scores_array = seq_trimmed_bio.quality_scores
    quality_scores_array.each do |each_qual|
      ee += 10**(-each_qual.to_f/10)
    end
    ee_2 = ee.round(2)
  end
    
  #puts ee
  return seq_trimmed.length, seq_trimmed, qual_trimmed, ee_2
end

##### Method to trim and orient the sequences
def orient (f_primer_matches, r_primer_matches, read_orientation, seq, qual)
  # Variable which will have the seq and qual stings oriented and trimmed
  seq_oriented = ""
  qual_oriented = ""

  # For cases when both primers matched
  if f_primer_matches == true and r_primer_matches == true
    if read_orientation == "+"
      #puts "#{f_primer_end}..#{r_primer_start}" 
      seq_oriented = seq
      qual_oriented = qual
    elsif read_orientation == "-"
      #puts "#{r_primer_end}..#{f_primer_start}"
      seq = Bio::Sequence::NA.new(seq)
      revcom_seq = seq.complement.upcase
      rev_qual = qual.reverse
      seq_oriented = revcom_seq
      qual_oriented = rev_qual
    end
 	end

  return seq_oriented, qual_oriented
end

##### Work with the file which has all the reads
def process_all_bc_reads_file (script_directory, all_bc_reads_file, ee, trim_req, human_db = nil, primer_file, thread, report_hash)
  # Get the basename of the fastq file
  file_basename = File.basename(all_bc_reads_file, ".*")

  # Opening the file in which all the information is going to be written to... and printing the header to it! 
  all_info_out_file = File.open("pre_all_reads_info.tsv", "w")
  all_info_out_file.puts("read_name\tbasename\tccs\tbarcode\tsample\tee_pretrim\tee_posttrim\tlength_pretrim\tlength_posttrim\thost_map\tf_primer_matches\tr_primer_matches\tf_primer_start\tf_primer_end\tr_primer_start\tr_primer_end\tread_orientation\tprimer_note\tnum_of_primer_hits")

  # Opening the file which which will have the trimmed and oriented sequences
  trimmed_out_file = File.open("pre_trimmed_and_oriented.fq", "w")

  # Get the seqs which map to the host genome by calling the map_to_host_genome method
  if !human_db.nil?
    puts "Mapping reads to reference genome...".green.bold
    mapped_count, mapped_string = map_to_human_genome(file_basename, human_db, thread) 
    puts "Done.".green.bold
  end

  # Primer matching by calling the primer_match method
  count_no_primer_match = 0
  puts "Aligning primers...".green.bold
  primer_record_hash = primer_match(script_directory, file_basename, primer_file, thread)
  abort("No primers found") if primer_record_hash.nil? || primer_record_hash.empty?
  puts "Done.".green.bold

  # Prereqs for the next loop
  #Create the hash which is going to store all info for each read using the Read_sequence class
  all_reads_hash = {}

  # Create a sequence hash which has all the sequence and quality strings for each record
  seqs_hash = {}

  # Create a hash with all the reads which are singletons
  singletons_hash = {}

  # Create a hash which has the trimmed seqs
  trimmed_hash = {}

  # Opening the file with all original reads for writing with bio module
  all_bc_reads = Bio::FlatFile.auto(all_bc_reads_file)
  
  puts "Parsing fastq headers and adding info to all_read table...".green.bold
  # Split out all attributes from the fastq header, check for required values 'ccs' and 'barcodelabel'
  all_bc_reads.each do |entry|
    def_split = entry.definition.split(";")
    all_rec_attr = {}
    read_name = def_split[0]
    def_split.each do |rec|
      if $. == 1
        all_rec_attr[:read_name] = rec 
      else
        var = rec.split("=")
        all_rec_attr[var[0]] = var[1]
      end
    end

    abort("Header does not have barcode label and ccs counts! Headers should have the format @read_name;barcodelabel=sample_name;ccs=ccs_passes. Can use get_fastqs.rb for this purpose.") unless all_rec_attr.has_key?("barcodelabel") && all_rec_attr.has_key?("ccs")    

    all_reads_hash[read_name]                = Read_sequence.new
    all_reads_hash[read_name].read_name      = all_rec_attr[:read_name]
    all_reads_hash[read_name].basename       = all_rec_attr["barcodelabel"]
    all_reads_hash[read_name].ccs            = all_rec_attr["ccs"].to_i
    all_reads_hash[read_name].barcode        = all_rec_attr["barcodelabel"]
    all_reads_hash[read_name].sample         = all_rec_attr["barcodelabel"]
    all_reads_hash[read_name].length_pretrim = entry.naseq.size

    #populating the report hash
    if report_hash.has_key?(all_reads_hash[read_name].barcode)
      report_hash[all_reads_hash[read_name].barcode].total += 1
    else 
      report_hash[all_reads_hash[read_name].barcode] = Report.new
      report_hash[all_reads_hash[read_name].barcode].total = 1
    end

    # Populate the seqs hash
    seqs_hash[read_name] = [entry.naseq.upcase, entry.quality_string, entry.definition]
  
    # Store the host genome mapping info
    if mapped_string.nil?
      all_reads_hash[read_name].host_map = "NA"
    elsif mapped_string.include?(read_name)
      all_reads_hash[read_name].host_map = true
    end 
    
    # Store the primer matching info
    if primer_record_hash.key?(read_name)
      if primer_record_hash[read_name][0] == true and primer_record_hash[read_name][1] == true and primer_record_hash[read_name][6] == "+" and primer_record_hash[read_name][3].to_i <= 100 and primer_record_hash[read_name][4].to_i >= all_reads_hash[read_name].length_pretrim-100
        all_reads_hash[read_name].f_primer_matches = primer_record_hash[read_name][0]
        all_reads_hash[read_name].r_primer_matches = primer_record_hash[read_name][1]
        all_reads_hash[read_name].f_primer_start = primer_record_hash[read_name][2].to_i
        all_reads_hash[read_name].f_primer_end = primer_record_hash[read_name][3].to_i
        all_reads_hash[read_name].r_primer_start = primer_record_hash[read_name][4].to_i
        all_reads_hash[read_name].r_primer_end = primer_record_hash[read_name][5].to_i
        all_reads_hash[read_name].read_orientation = primer_record_hash[read_name][6]
        all_reads_hash[read_name].primer_note = primer_record_hash[read_name][8]
        all_reads_hash[read_name].num_of_primerhits = primer_record_hash[read_name][9]
      elsif primer_record_hash[read_name][0] == true and primer_record_hash[read_name][1] == true and primer_record_hash[read_name][6] == "-" and primer_record_hash[read_name][2].to_i >= all_reads_hash[read_name].length_pretrim-100 and primer_record_hash[read_name][5].to_i <= 100
        all_reads_hash[read_name].f_primer_matches = primer_record_hash[read_name][0]
        all_reads_hash[read_name].r_primer_matches = primer_record_hash[read_name][1]
        all_reads_hash[read_name].f_primer_start = primer_record_hash[read_name][2].to_i
        all_reads_hash[read_name].f_primer_end = primer_record_hash[read_name][3].to_i
        all_reads_hash[read_name].r_primer_start = primer_record_hash[read_name][4].to_i
        all_reads_hash[read_name].r_primer_end = primer_record_hash[read_name][5].to_i
        all_reads_hash[read_name].read_orientation = primer_record_hash[read_name][6]
        all_reads_hash[read_name].primer_note = primer_record_hash[read_name][8]
        all_reads_hash[read_name].num_of_primerhits = primer_record_hash[read_name][9]
      elsif primer_record_hash[read_name][0] == false or primer_record_hash[read_name][1] == false
        all_reads_hash[read_name].primer_note = primer_record_hash[read_name][8]
        all_reads_hash[read_name].f_primer_matches = primer_record_hash[read_name][0]
        all_reads_hash[read_name].r_primer_matches = primer_record_hash[read_name][1]
        all_reads_hash[read_name].f_primer_start = primer_record_hash[read_name][2].to_i
        all_reads_hash[read_name].f_primer_end = primer_record_hash[read_name][3].to_i
        all_reads_hash[read_name].r_primer_start = primer_record_hash[read_name][4].to_i
        all_reads_hash[read_name].r_primer_end = primer_record_hash[read_name][5].to_i
        all_reads_hash[read_name].read_orientation = primer_record_hash[read_name][6]
        all_reads_hash[read_name].num_of_primerhits = primer_record_hash[read_name][9]
       	singletons_hash[read_name] = [primer_record_hash[read_name][0], primer_record_hash[read_name][1]]
      else
        all_reads_hash[read_name].primer_note = "primer_match_length_criteria_failed"
      end
    else
      count_no_primer_match += 1
      all_reads_hash[read_name].primer_note = "no_primer_hits"
    end
  end

  puts "Done.".green.bold
  
  # Get ee_pretrim
  puts "Calculating EE for all fastq sequences...".green.bold
  ee_pretrim_hash = get_ee_from_fq_file(file_basename, ee, "pre_ee_pretrim.fq")
  puts "Done.".green.bold

  ee_pretrim_hash.each do |k, v|
    all_reads_hash[k].ee_pretrim = v
  end

  # Add all remaining values to to each read record (trimming info, ee, length etc.)
  puts "Parsing ee and trimming info for every read...".green.bold
  all_reads_hash.each do |k, v|
    if trim_req == "yes"
    	# Trim and orient sequences with the trim_and_orient method
    	seq_trimmed_length, seq_trimmed, qual_trimmed, ee = trim_and_orient(all_reads_hash[k].f_primer_matches, all_reads_hash[k].r_primer_matches, all_reads_hash[k].f_primer_start, all_reads_hash[k].f_primer_end, all_reads_hash[k].r_primer_start, all_reads_hash[k].r_primer_end, all_reads_hash[k].read_orientation, seqs_hash[k][0], seqs_hash[k][1])
    	if seq_trimmed == ""
      	# Add the length_postrim info to the all_reads_hash
      	all_reads_hash[k].length_posttrim = -1
      	# Add the ee postrim to the all_reads_hash
      	all_reads_hash[k].ee_posttrim = -1
    	else
      	# Add the length_postrim info to the all_reads_hash
      	all_reads_hash[k].length_posttrim = seq_trimmed_length
      	# Add the ee postrim to the all_reads_hash
      	all_reads_hash[k].ee_posttrim = ee
      	# write the oriented and trimmed sequences to an output file
      	write_to_fastq(trimmed_out_file, seqs_hash[k][2], seq_trimmed, qual_trimmed)
      	# populate the hash havng the trimmed seqs
      	trimmed_hash[seqs_hash[k][2]] = [seq_trimmed, qual_trimmed]
    	end
    else
    	# No trimming, only orient sequences using the orient method
			seq_oriented, qual_oriented = orient(all_reads_hash[k].f_primer_matches, all_reads_hash[k].r_primer_matches, all_reads_hash[k].read_orientation, seqs_hash[k][0], seqs_hash[k][1])
    	if seq_oriented == ""
      	# Add the length_postrim info to the all_reads_hash
      	all_reads_hash[k].length_posttrim = -1
      	# Add the ee postrim to the all_reads_hash
      	all_reads_hash[k].ee_posttrim = -1
    	else
      	# Add the length_postrim info to the all_reads_hash
      	all_reads_hash[k].length_posttrim = seq_oriented.length
      	# Add the ee postrim to the all_reads_hash
        all_reads_hash[k].ee_posttrim = all_reads_hash[k].ee_pretrim
      	# write the oriented and trimmed sequences to an output file
      	write_to_fastq(trimmed_out_file, seqs_hash[k][2], seq_oriented, qual_oriented)
      	# populate the hash havng the trimmed seqs
      	trimmed_hash[seqs_hash[k][2]] = [seq_oriented, qual_oriented]
    	end
    end

    # Write all reads to output file
    all_info_out_file.puts([k,
                            all_reads_hash[k].basename,
                            all_reads_hash[k].ccs,
                            all_reads_hash[k].barcode,
                            all_reads_hash[k].sample,
                            all_reads_hash[k].ee_pretrim,
                            all_reads_hash[k].ee_posttrim,
                            all_reads_hash[k].length_pretrim,
                            all_reads_hash[k].length_posttrim,
                            all_reads_hash[k].host_map,
                            all_reads_hash[k].f_primer_matches,
                            all_reads_hash[k].r_primer_matches,
                            all_reads_hash[k].f_primer_start,
                            all_reads_hash[k].f_primer_end,
                            all_reads_hash[k].r_primer_start,
                            all_reads_hash[k].r_primer_end,
                            all_reads_hash[k].read_orientation,
                            all_reads_hash[k].primer_note,
                            all_reads_hash[k].num_of_primerhits].join("\t"))		
  end
  puts "Done.".green.bold
  # Close the output file in which the trimmed and oriented seqs were written
  trimmed_out_file.close

  # Close the output file which has all the info
  all_info_out_file.close   

  # return the all info hash and trimmed hash
  return all_reads_hash, trimmed_hash, report_hash
end


##################### MAIN PROGRAM #######################

# Calling the method which comcatenates files
puts "Concatenating fastq files..".green.bold
concat_files(folder_name, all_files, sample_list)
puts "Done.".green.bold

# Getting the name of the file which has all the reads
all_bc_reads_file = "pre_demultiplexed_ccsfilt.fq"

# Create a hash which populates the Report class
report_hash = {}

# Create a file which will have the report for number of sequences in each step
puts "Initializing report file...".green.bold
report_file = File.open("post_each_step_report.tsv", "w") 
# writing into the report file
report_file.puts("sample\tdemultiplexed_ccs\tprimer_filt\thost_filt\tsize_filt\tee_filt\tderep_filt\tchimera_filt")
puts "Done.".green.bold

# Calling the method which then calls all the other methods! 
puts "Starting main method...".magenta.bold
all_reads_hash, trimmed_hash, report_hash = process_all_bc_reads_file(script_directory, all_bc_reads_file, ee, trim_req, human_db, primer_file, thread, report_hash)
puts "Finished main method.".magenta.bold

final_fastq_file = File.open("pre_sequences_for_clustering.fq", "w")
final_fastq_basename = File.basename(final_fastq_file, ".fq")

file_for_usearch_global = File.open("pre_sequences_for_usearch_global.fq", "w")
file_for_usearch_global_basename = File.basename(file_for_usearch_global)

# loop through the hash which has the trimmed sequences
puts "Writing trimmed fastq values to file...".green.bold
trimmed_hash.each do |key, value|
  #puts key, value
  key_in_all_reads_hash = key.split(";")[0]
  #puts all_reads_hash[key_in_all_reads_hash].read_name

  # Writing to files based on filtering
  if all_reads_hash[key_in_all_reads_hash].ccs >= ccs and all_reads_hash[key_in_all_reads_hash].ee_posttrim <= ee and all_reads_hash[key_in_all_reads_hash].length_posttrim >= length_min and all_reads_hash[key_in_all_reads_hash].length_posttrim <= length_max and all_reads_hash[key_in_all_reads_hash].host_map != true
    write_to_fastq(final_fastq_file, key, value[0], value[1])
  end

  if all_reads_hash[key_in_all_reads_hash].ccs >= ccs and all_reads_hash[key_in_all_reads_hash].length_posttrim >= length_min and all_reads_hash[key_in_all_reads_hash].length_posttrim <= length_max and all_reads_hash[key_in_all_reads_hash].host_map != true
    write_to_fastq(file_for_usearch_global, key, value[0], value[1])
  end

  # Populating the report_hash using the report class
  basename_key = all_reads_hash[key_in_all_reads_hash].basename
  if all_reads_hash[key_in_all_reads_hash].primer_note == "good"
    report_hash[basename_key].primer_filt += 1
  end
  if all_reads_hash[key_in_all_reads_hash].primer_note == "good" and all_reads_hash[key_in_all_reads_hash].host_map != true
    report_hash[basename_key].host_filt += 1
  end
  if all_reads_hash[key_in_all_reads_hash].primer_note == "good" and all_reads_hash[key_in_all_reads_hash].host_map != true and all_reads_hash[key_in_all_reads_hash].length_posttrim >= length_min and all_reads_hash[key_in_all_reads_hash].length_posttrim <= length_max
    report_hash[basename_key].size_filt += 1
  end
  if all_reads_hash[key_in_all_reads_hash].primer_note == "good" and all_reads_hash[key_in_all_reads_hash].host_map != true and all_reads_hash[key_in_all_reads_hash].length_posttrim >= length_min and all_reads_hash[key_in_all_reads_hash].length_posttrim <= length_max and all_reads_hash[key_in_all_reads_hash].ee_posttrim <= ee
    report_hash[basename_key].ee_filt += 1
  end

end
puts "Done.".green.bold

final_fastq_file.close
file_for_usearch_global.close

# Running the usearch commands for clustering
puts "Clustering...".magenta.bold
`sh #{script_directory}/uparse_commands.sh #{final_fastq_basename} #{uchime_db_file} #{utax_db_file} #{file_for_usearch_global_basename}`
puts "Done".magenta.bold

# Running the command to give a report of counts
puts "Generating Reports...".green.bold
`ruby #{script_directory}/get_report.rb #{final_fastq_basename}`
puts "Done.".green.bold

# Running blast on the OTUs
if !lineage_fasta_file.nil?
  puts "Blasting OTU centroids...".magenta.bold                                                                                  
  `usearch -ublast post_OTU.fa -db #{lineage_fasta_file} -top_hit_only -id 0.9 -blast6out post_OTU_blast.txt -strand both -evalue 0.01 -threads #{thread} -accel 0.3`
  puts "Done.".magenta.bold

  # Running the script whcih gives a final file with all the clustering info, taxa info and blast info
  `ruby #{script_directory}/final_parsing.rb -b post_OTU_blast.txt -u post_OTU_table_utax_map.tsv -n #{ncbi_clust_file} -o post_final_results.tsv`
else  
  `ruby #{script_directory}/final_parsing.rb -u post_OTU_table_utax_map.tsv -n #{ncbi_clust_file} -o post_final_results.tsv`
end


# Run usearch on all the reads
puts "Utaxing cluster centroids...".green.bold
`usearch -utax #{all_bc_reads_file} -db #{utax_db_file} -utaxout all_bc_reads.utax -utax_cutoff 0.8 -strand both -threads #{thread}`
puts "Done.".green.bold

# Run ruby script to merge the all_bc_info file and all_bc_utax files
`ruby #{script_directory}/merge_all_info_and_utax.rb`    

if(File.directory?('split_otus')) 
  `rm -rf split_otus`
  `mkdir split_otus`
else
  `mkdir split_otus` 
end

if split_otus == "yes" and split_otu_method == "before"
  `ruby #{script_directory}/extract_seqs_for_each_otu_noEEfilt.rb -u post_OTU_nonchimeras.fa -p post_usearch_glob_results.tsv -a pre_sequences_for_usearch_global.fq`
elsif split_otus == "yes" and split_otu_method == "after"
  `ruby #{script_directory}/extract_seqs_for_each_otu.rb -u post_OTU_nonchimeras.fa -p post_uparse.up -a post_dereplicated.fa`
end 

# Dealing with the dereplicated file 
derep_open =  Bio::FlatFile.auto("post_dereplicated.fa")                                                                                                     
derep_open.each do |entry|
  def_split = entry.definition.split(";")
  basename_ind = def_split.index{|s| s.include?("barcodelabel")}
  basename = def_split[basename_ind].split("=")[1]
  report_hash[basename].derep_filt += 1
end

# Dealing with the up file from clustering step
uparse_open = File.open("post_uparse.up")
uparse_open.each do |line|
  line_split = line.split("\t")[0].split(";")
  basename_ind = line_split.index{|s| s.include?("barcodelabel")}
  basename = line_split[basename_ind].split("=")[1]
  #puts basename
  type = line_split[1]
  if type != "chimera"
    report_hash[basename].chimera_filt += 1
  end
end

# and since that is the last thing to be added to the report hash, print it out

#puts report_hash
report_hash.each do |key, value|
  #puts key, value.total
  report_file.puts([key,
                    value.total, 
                    value.primer_filt,
                    value.host_filt,
                    value.size_filt,
                    value.ee_filt,
                    value.derep_filt,
                    value.chimera_filt].join("\t"))
end

if verbose == true
  abort
else
  File.delete("post_dereplicated.fa")
  File.delete("post_OTU_alignment.aln")
  File.delete("post_OTU_blast.txt")
  File.delete("post_OTU_chimeras.fa")
  File.delete("post_OTU_nonchimeras.fa")
  File.delete("post_OTU_table.tsv")
  File.delete("post_OTU_table_utax_map.tsv")
  File.delete("post_OTU_uchime_output.tsv")
  File.delete("post_readmap.uc")
  File.delete("post_reads.utax")
  File.delete("post_unmapped_userach_global.fa")
  File.delete("post_uparse.up")
  File.delete("pre_ee_pretrim.fq")
  File.delete("pre_filt_non_host.fq")
  File.delete("pre_host_mapped.txt")
  File.delete("pre_map_to_host.bam")
  File.delete("pre_map_to_host.sam")
  File.delete("pre_primer_map_info.tsv")
  File.delete("pre_primer_map.tsv")
  File.delete("pre_sequences_for_clustering.fq")
  File.delete("pre_sequences_for_usearch_global.fq")
  File.delete("pre_trimmed_and_oriented.fq")
end
#=end
