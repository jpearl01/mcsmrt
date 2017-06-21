# -*- coding: utf-8 -*-
require 'bio'
require 'trollop'

#USAGE: ruby ../mcsmrt/mcsmrt_v1.rb -a -f reads/ -d 32 -e 1 -s 5 -x 2000 -n 500 -c ../rdp_gold.fa -t ../lineanator/16sMicrobial_ncbi_lineage_reference_database.udb -l ../lineanator/16sMicrobial_ncbi_lineage.fasta -g ../human_g1k_v37.fasta -p ../primers.fasta

##### Input 
opts = Trollop::options do
  opt :allFiles, "Use all files in the given directory name for clustering", :short => "-a"
  opt :foldername, "Folder with the demultiplexed files for clustering", :type => :string, :short => "-f"
  opt :samplelist, "File with a list of file names which are to be merged and clustered together", :type => :string, :short => "-i"
  opt :threads, "Number of threads you can allot for running this process", :type => :int, :short => "-d", :default => 1
  opt :eevalue, "Expected error value greater than which reads will be filtered out", :type => :float, :short => "-e", :default => 1.0
  opt :trimming, "Do you want to trim your sequences? Answer in yes or no", :type => :string, :short => "-m", :default => "yes"
  opt :ccsvalue, "CCS passes lesser than which reads will be filtered out", :type => :int, :short => "-s", :default => 5
  opt :lengthmax, "Maximum length above which reads will be filtered out", :type => :int, :short => "-x", :default => 2000
  opt :lengthmin, "Manimum length below which reads will be filtered out", :type => :int, :short => "-n", :default => 500
  opt :uchimedbfile, "Path to database file for the uchime command", :type => :string, :short => "-c"
  opt :utaxdbfile, "Path to database file for the utax command", :type => :string, :short => "-t"
  opt :lineagefastafile, "Path to FASTA file with lineage info for the ublast command", :type => :string, :short => "-l"
  opt :host_db, "Path to fasta file of host genome", :type => :string, :short => "-g"
  opt :primerfile, "Path to fasta file with the primer sequences", :type => :string, :short => "-p"
end 

##### Assigning variables to the input and making sure we got all the inputs
if opts[:allFiles] == true     
  all_files = "*.fq"  
else
  all_files = nil   
end  

if opts[:samplelist] != nil  
  sample_list = opts[:samplelist]
else
  sample_list = nil
end              

if opts[:allFiles] == false and opts[:samplelist] == nil
  abort("Must specify if you want to use all files in the folder with '-a' or give a file with a list of file names for clustering with '-i'")
end

opts[:foldername].nil?        ==false  ? folder_name = opts[:foldername]              : abort("Must supply the name of the folder in which the demultiplexed files exist with '-f'")
opts[:uchimedbfile].nil?      ==false  ? uchime_db_file = opts[:uchimedbfile]         : abort("Must supply a 'uchime database file' e.g. rdpgold.udb with '-c'")
opts[:utaxdbfile].nil?        ==false  ? utax_db_file = opts[:utaxdbfile]             : abort("Must supply a 'utax database file' e.g. 16s_ncbi.udb with '-t'")
opts[:lineagefastafile].nil?  ==false  ? lineage_fasta_file = opts[:lineagefastafile] : abort("Must supply a 'lineage fasta file' e.g. ncbi_lineage.fasta (for blast) with '-l'")
opts[:host_db].nil?           ==false  ? human_db = opts[:host_db]                    : abort("Must supply a fasta of the host genome e.g. human_g1k.fasta with '-g'")
opts[:primerfile].nil?        ==false  ? primer_file = opts[:primerfile]              : abort("Must supply a fasta of the primer sequences e.g primer_seqs.fa with '-p'")

thread = opts[:threads].to_i
ee = opts[:eevalue].to_f 
trim_req = opts[:trimming].to_s
ccs = opts[:ccsvalue].to_i 
length_max = opts[:lengthmax].to_i 
length_min = opts[:lengthmin].to_i 

##### Get the path to the directory in which the scripts exist 
script_directory = File.dirname(__FILE__)

##### Class that stores information about each record from the reads file
class Read_sequence
  attr_accessor :read_name, :basename, :ccs, :barcode, :sample, :ee_pretrim, :ee_posttrim, :length_pretrim, :length_posttrim, :host_map, :f_primer_matches, :r_primer_matches, :f_primer_start, :f_primer_end, :r_primer_start, :r_primer_end, :read_orientation, :primer_note

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
    @host_map          = false 
    @f_primer_matches  = false
    @r_primer_matches  = false
    @f_primer_start    = "NA"
    @f_primer_end      = "NA"
    @r_primer_start    = "NA"
    @r_primer_end      = "NA"
    @read_orientation  = "NA"
    @primer_note       = "NA"
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
    `cat #{list.join(" ")} > demultiplexed_ccsfilt.fq`

  else
    `cat #{folder_name}/#{all_files} > demultiplexed_ccsfilt.fq` 
  end

  abort("!!!!The file with all the reads required for clustering does not exist!!!!") if !File.exists?("demultiplexed_ccsfilt.fq")
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
	`usearch -fastq_filter #{file_basename}.fq -fastqout #{file_basename}_#{suffix} -fastq_maxee 20000 -fastq_qmax 75 -fastq_eeout -sample all`
	
  abort("!!!!Expected error filtering with usearch failed!!!!") if File.zero?("#{file_basename}_#{suffix}")

	# Hash that is returned from this method (read name - key, ee - value)
	ee_hash = {}
	
	# Open the file from fastq_filter command and process it
	ee_filt_file = Bio::FlatFile.auto("#{file_basename}_#{suffix}")
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
  `bwa mem -t #{thread} #{human_db} #{file_basename}.fq > map_to_host.sam`

  # Check to make sure bwa worked, which means that the index files also exist
  abort("!!!!Index files are missing, run the bwa index command to index the reference genome and store them in the same directory as the reference genome!!!!") if File.zero?("map_to_host.sam")
  
  #sambamba converts sam to bam format                                                                                                                   
  `sambamba view -S -f bam map_to_host.sam -o map_to_host.bam`
  
  #Sort the bam file                                                                                                                                     
  #`sambamba sort -t#{thread} -o filt_non_host_sorted.bam filt_non_host.bam`
  
  #filter the bam for only ‘not unmapped’ reads -> reads that are mapped                                                                                 
  `sambamba view -F 'not unmapped' map_to_host.bam > host_mapped.txt`
  mapped_count = `cut -d ';' -f1 host_mapped.txt | sort | uniq | wc -l`
  mapped_string = `cut -d ';' -f1 host_mapped.txt`
  
  #filter reads out for ‘unmapped’ -> we would use these for pipeline                                                                             
  `sambamba view -F 'unmapped' map_to_host.bam > filt_non_host.txt`
  
  #convert the sam file to fastq                                                                                                                         
  `grep -v ^@ filt_non_host.txt | awk '{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}' > filt_non_host.fq`
  
 	return mapped_count, mapped_string
end

##### Method for primer matching 
def primer_match (script_directory, file_basename, primer_file, thread)
	# Run the usearch command for primer matching
	`usearch -usearch_local #{file_basename}.fq -db #{primer_file} -id 0.8 -threads #{thread} -userout primer_map.txt -userfields query+target+qstrand+tlo+thi+qlo+qhi+pairs+gaps+mism+diffs -strand both -gapopen 1.0 -gapext 0.5 -lopen 1.0 -lext 0.5`
  #`usearch -search_oligodb #{file_basename}.fq -db #{primer_file} -strand both -userout #{file_basename}_primer_map.txt -userfields query+target+qstrand+diffs+tlo+thi+qlo+qhi`                                                                                                    

  # Check to see if sequences passed primer matching, i.e., if no read has a hit for primers, the output file from the previous step will be empty!
  abort("!!!!None of the reads mapped to the primers, check your FASTA file which has the primers!!!!") if File.zero?("primer_map.txt")

  # Run the script which parses the primer matching output
  `ruby #{script_directory}/primer_matching_usearch_local.rb -p primer_map.txt -o primer_map_info.txt -a #{file_basename}.fq` 
  	
  # Open the file with parsed primer matching results
  "primer_map_info.txt".nil? ==false  ? primer_matching_parsed_file = File.open("primer_map_info.txt") : abort("Primer mapping parsed file was not created from primer_matching.rb")

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
      #puts "#{f_primer_end}..#{r_primer_start}" 
      seq_trimmed = seq[f_primer_end..r_primer_start]
      qual_trimmed = qual[f_primer_end..r_primer_start]
    elsif read_orientation == "-"
      #puts "#{r_primer_end}..#{f_primer_start}"
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
      prob_wrong = -each_qual.to_f/10
      prob_wrong_2 = 10**prob_wrong
      ee += prob_wrong_2
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
def process_all_bc_reads_file (script_directory, all_bc_reads_file, ee, trim_req, human_db, primer_file, thread)
  # Get the basename of the fastq file
  file_basename = File.basename(all_bc_reads_file, ".*")

  # Opening the file in which all the information is going to be written to... and printing the header to it! 
  all_info_out_file = File.open("all_reads_info.txt", "w")
  all_info_out_file.puts("read_name\tbasename\tccs\tbarcode\tsample\tee_pretrim\tee_posttrim\tlength_pretrim\tlength_posttrim\thost_map\tf_primer_matches\tr_primer_matches\tf_primer_start\tf_primer_end\tr_primer_start\tr_primer_end\tread_orientation\tprimer_note")

  # Opening the file which which will have the trimmed and oriented sequences
  trimmed_out_file = File.open("trimmed_and_oriented.fq", "w")

  # Get the seqs which map to the host genome by calling the map_to_host_genome method
 	mapped_count, mapped_string = map_to_human_genome(file_basename, human_db, thread) 
  #puts mapped_string.inspect 

  # Primer matching by calling the primer_match method
  count_no_primer_match = 0
  primer_record_hash = primer_match(script_directory, file_basename, primer_file, thread)
  #puts primer_record_hash.inspect
  #puts primer_record_hash["m150212_113119_42168_c100747522550000001823168507081543_s1_p0/32770/ccs"]
  abort

  # Prereqs for the next loop
 	#Create the hash which is going to store all infor for each read using the Read_sequence class
  all_reads_hash = {}
  # Create a sequence hash which has all the sequence and quality strings for each record
  seqs_hash = {}
  # Create a hash with all the reads which are singletons
  singletons_hash = {}
  # Create a hash which has the trimmed seqs
	trimmed_hash = {}
  # Opening the file with all original reads for writing with bio module
  all_bc_reads = Bio::FlatFile.auto(all_bc_reads_file)

  # Loop through the fq file, record some basic infor in the all_reads_hash
  all_bc_reads.each do |entry|
    def_split = entry.definition.split(";")
    read_name = def_split[0]
    #puts def_split

    # Fill the all_bc_hash with some basic info that we can get (read_name, barcode, ccs, length_pretrim)
    if def_split[1].include?("barcodelabel")
      basename = def_split[1].split("=")[1]
      barcode = basename.split("_")[0]
      sample = basename.split("_")[1..-1].join("_")
      ccs = def_split[2].split("=")[1].to_i
    elsif def_split[1].include?("ccs")
      basename = def_split[2].split("=")[1]
      barcode = basename.split("_")[0]
      sample = basename.split("_")[1..-1].join("_")
      ccs = def_split[1].split("=")[1].to_i
    else
      abort("Header does not have barcode label and ccs counts! Headers should have the format @read_name;barcodelabel=sample_name;ccs=ccs_passes. Use get_fastqs.rb for this purpose.")
    end

    all_reads_hash[read_name] = Read_sequence.new
    all_reads_hash[read_name].read_name = read_name
    all_reads_hash[read_name].basename = basename
    all_reads_hash[read_name].ccs = ccs
    all_reads_hash[read_name].barcode = barcode
    all_reads_hash[read_name].sample = sample
    all_reads_hash[read_name].length_pretrim = entry.naseq.size

  	# Populate the seqs hash
    seqs_hash[read_name] = [entry.naseq.upcase, entry.quality_string, entry.definition]
  
    # Store the host genome mapping info
    if mapped_string.include?(read_name)
      all_reads_hash[read_name].host_map = true
    end 
    
    # Store the primer matching info
    if primer_record_hash.key?(read_name)
      #puts primer_record_hash[read_name].inspect
      #puts primer_record_hash[read_name][0], primer_record_hash[read_name][1], primer_record_hash[read_name][6], primer_record_hash[read_name][3].to_i, primer_record_hash[read_name][4].to_i, primer_record_hash[read_name][2].to_i, primer_record_hash[read_name][5].to_i
      #puts all_reads_hash[read_name].length_pretrim.inspect
      if primer_record_hash[read_name][0] == true and primer_record_hash[read_name][1] == true and primer_record_hash[read_name][6] == "+" and primer_record_hash[read_name][3].to_i <= 100 and primer_record_hash[read_name][4].to_i >= all_reads_hash[read_name].length_pretrim-100
        all_reads_hash[read_name].f_primer_matches = primer_record_hash[read_name][0]
        all_reads_hash[read_name].r_primer_matches = primer_record_hash[read_name][1]
        all_reads_hash[read_name].f_primer_start = primer_record_hash[read_name][2].to_i
        all_reads_hash[read_name].f_primer_end = primer_record_hash[read_name][3].to_i
        all_reads_hash[read_name].r_primer_start = primer_record_hash[read_name][4].to_i
        all_reads_hash[read_name].r_primer_end = primer_record_hash[read_name][5].to_i
        all_reads_hash[read_name].read_orientation = primer_record_hash[read_name][6]
        all_reads_hash[read_name].primer_note = primer_record_hash[read_name][8]
      elsif primer_record_hash[read_name][0] == true and primer_record_hash[read_name][1] == true and primer_record_hash[read_name][6] == "-" and primer_record_hash[read_name][2].to_i >= all_reads_hash[read_name].length_pretrim-100 and primer_record_hash[read_name][5].to_i <= 100
        all_reads_hash[read_name].f_primer_matches = primer_record_hash[read_name][0]
        all_reads_hash[read_name].r_primer_matches = primer_record_hash[read_name][1]
        all_reads_hash[read_name].f_primer_start = primer_record_hash[read_name][2].to_i
        all_reads_hash[read_name].f_primer_end = primer_record_hash[read_name][3].to_i
        all_reads_hash[read_name].r_primer_start = primer_record_hash[read_name][4].to_i
        all_reads_hash[read_name].r_primer_end = primer_record_hash[read_name][5].to_i
        all_reads_hash[read_name].read_orientation = primer_record_hash[read_name][6]
        all_reads_hash[read_name].primer_note = primer_record_hash[read_name][8]
      elsif primer_record_hash[read_name][0] == false or primer_record_hash[read_name][1] == false
        all_reads_hash[read_name].primer_note = primer_record_hash[read_name][8]
        all_reads_hash[read_name].f_primer_matches = primer_record_hash[read_name][0]
        all_reads_hash[read_name].r_primer_matches = primer_record_hash[read_name][1]
        all_reads_hash[read_name].f_primer_start = primer_record_hash[read_name][2].to_i
        all_reads_hash[read_name].f_primer_end = primer_record_hash[read_name][3].to_i
        all_reads_hash[read_name].r_primer_start = primer_record_hash[read_name][4].to_i
        all_reads_hash[read_name].r_primer_end = primer_record_hash[read_name][5].to_i
        all_reads_hash[read_name].read_orientation = primer_record_hash[read_name][6]
       	singletons_hash[read_name] = [primer_record_hash[read_name][0], primer_record_hash[read_name][1]]
      else
        all_reads_hash[read_name].primer_note = "primer_match_length_criteria_failed"
      end
    else
      count_no_primer_match += 1
      all_reads_hash[read_name].primer_note = "no_primer_hits"
    end
  end
  
  # Get ee_pretrim
  ee_pretrim_hash = get_ee_from_fq_file(file_basename, ee, "ee_pretrim.fq")
  #puts ee_pretrim_hash
  ee_pretrim_hash.each do |k, v|
   	#puts all_reads_hash[read_name]
    all_reads_hash[k].ee_pretrim = v
 	end
  
  #puts all_reads_hash[k].f_primer_matches, all_reads_hash[k].r_primer_matches, all_reads_hash[k].f_primer_start, all_reads_hash[k].f_primer_end, all_reads_hash[k].r_primer_start, all_reads_hash[k].r_primer_end, all_reads_hash[k].read_orientation, all_reads_hash[k].half_primer_match, all_reads_hash[k].half_match_start, all_reads_hash[k].half_match_end, seqs_hash[k][0], seqs_hash[k][1]
	
  # Loop thorugh the all_reads_hash and fill the rest of the info into it
	all_reads_hash.each do |k, v|
    if trim_req == "yes"
    	# Trim and orient sequences with the trim_and_orient method
    	seq_trimmed_length, seq_trimmed, qual_trimmed, ee = trim_and_orient(all_reads_hash[k].f_primer_matches, all_reads_hash[k].r_primer_matches, all_reads_hash[k].f_primer_start, all_reads_hash[k].f_primer_end, all_reads_hash[k].r_primer_start, all_reads_hash[k].r_primer_end, all_reads_hash[k].read_orientation, seqs_hash[k][0], seqs_hash[k][1])
    	#puts "#{seq_trimmed_length},#{seq_trimmed},#{qual_trimmed}" 
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
    all_info_out_file.puts("#{k}\t#{all_reads_hash[k].basename}\t#{all_reads_hash[k].ccs}\t#{all_reads_hash[k].barcode}\t#{all_reads_hash[k].sample}\t#{all_reads_hash[k].ee_pretrim}\t#{all_reads_hash[k].ee_posttrim}\t#{all_reads_hash[k].length_pretrim}\t#{all_reads_hash[k].length_posttrim}\t#{all_reads_hash[k].host_map}\t#{all_reads_hash[k].f_primer_matches}\t#{all_reads_hash[k].r_primer_matches}\t#{all_reads_hash[k].f_primer_start}\t#{all_reads_hash[k].f_primer_end}\t#{all_reads_hash[k].r_primer_start}\t#{all_reads_hash[k].r_primer_end}\t#{all_reads_hash[k].read_orientation}\t#{all_reads_hash[k].primer_note}")		
  end

  # Close the output file in which the trimmed and oriented seqs were written
  trimmed_out_file.close

  # Close the output file which has all the info
  all_info_out_file.close   

  # return the all info hash and trimmed hash
  return all_reads_hash, trimmed_hash 
end


##################### MAIN PROGRAM #######################

# Calling the method which comcatenates files
concat_files(folder_name, all_files, sample_list)

# Getting the name of the file which has all the reads
all_bc_reads_file = "demultiplexed_ccsfilt.fq"

# Calling the method which then calls all the other methods! 
all_reads_hash, trimmed_hash = process_all_bc_reads_file(script_directory, all_bc_reads_file, ee, trim_req, human_db, primer_file, thread)

#=begin

final_fastq_file = File.open("sequences_for_clustering.fq", "w")
final_fastq_basename = File.basename(final_fastq_file, ".fq")

file_for_usearch_global = File.open("sequences_for_usearch_global.fq", "w")
file_for_usearch_global_basename = File.basename(file_for_usearch_global)

trimmed_hash.each do |key, value|
  #puts key, value
  key_in_all_reads_hash = key.split(";")[0]
  #puts all_reads_hash[key_in_all_reads_hash].read_name

  if all_reads_hash[key_in_all_reads_hash].ccs >= ccs and all_reads_hash[key_in_all_reads_hash].ee_posttrim <= ee and all_reads_hash[key_in_all_reads_hash].length_posttrim >= length_min and all_reads_hash[key_in_all_reads_hash].length_posttrim <= length_max and all_reads_hash[key_in_all_reads_hash].host_map == false
    write_to_fastq(final_fastq_file, key, value[0], value[1])
    #else
    #  puts all_reads_hash[key_in_all_reads_hash].ccs, all_reads_hash[key_in_all_reads_hash].ee_posttrim, all_reads_hash[key_in_all_reads_hash].length_posttrim, all_reads_hash[key_in_all_reads_hash].host_map
    #end
  end

  if all_reads_hash[key_in_all_reads_hash].ccs >= ccs and all_reads_hash[key_in_all_reads_hash].length_posttrim >= length_min and all_reads_hash[key_in_all_reads_hash].length_posttrim <= length_max and all_reads_hash[key_in_all_reads_hash].host_map == false
    write_to_fastq(file_for_usearch_global, key, value[0], value[1])
  end

end

final_fastq_file.close
file_for_usearch_global.close

# Running the usearch commands for clustering
`sh #{script_directory}/uparse_commands.sh #{final_fastq_basename} #{uchime_db_file} #{utax_db_file} #{file_for_usearch_global_basename}`

# Running the command to give a report of counts
`ruby #{script_directory}/get_report.rb #{final_fastq_basename}`
  
# Running blast on the OTUs                                                                                                                            
`usearch -ublast #{final_fastq_basename}_OTU_s2.fa -db #{lineage_fasta_file} -top_hit_only -id 0.9 -blast6out #{final_fastq_basename}_blast.txt -strand both -evalue 0.01 -threads #{thread} -accel 0.3`

# Running the script whcih gives a final file with all the clustering info, taxa info and blast info
`ruby #{script_directory}/final_parsing.rb -b #{final_fastq_basename}_blast.txt -u #{final_fastq_basename}_OTU_table_utax_map.txt -o #{final_fastq_basename}_final.txt`

#=end
