#!/usr/bin/env ruby
# -*- coding: utf-8 -*-

require 'bio'
require 'trollop'


#USAGE: ruby pool_samples.rb -s sample_key_BEI.txt -e 1 -c ../rdp_gold.fa -t ../16s_ncbi_database/16s_lineage_short_species_name_reference_database.udb -l ../16s_ncbi_database/16sMicrobial_ncbi_lineage.fasta

##### Input 
opts = Trollop::options do
	opt :samplefile, "File with all the sample information", :type => :string, :short => "-s"
	opt :eevalue, "Expected error at which you want filtering to take place", :type => :string, :short => "-e"
	opt :uchimedbfile, "Path of database file for the uchime command", :type => :string, :short => "-c"
	opt :utaxdbfile, "Path of database file for the utax command", :type => :string, :short => "-t"
	opt :lineagefastafile, "Path of FASTA file with lineage info for the ublast command", :type => :string, :short => "-l"
end 

##### Assigning variables to the input
opts[:samplefile].nil?       ==false  ? sample_file = opts[:samplefile]               : abort("Must supply a 'sample file': tab delimited file of sample information with '-s'")
opts[:eevalue].nil?          ==false  ? ee = opts[:eevalue]                           : abort("Must supply an Expected Error value with '-e'")
opts[:uchimedbfile].nil?     ==false  ? uchime_db_file = opts[:uchimedbfile]          : abort("Must supply a 'uchime database file' e.g. rdpgold.udb '-c'")
opts[:utaxdbfile].nil?       ==false  ? utax_db_file = opts[:utaxdbfile]              : abort("Must supply a 'utax database file' e.g. 16s_ncbi.udb '-t'")
opts[:lineagefastafile].nil? ==false  ? lineage_fasta_file = opts[:lineagefastafile]  : abort("Must supply a 'lineage fasta file' e.g. ncbi_lineage.fasta (for blast) '-l' ")

##### Making sure we got the inputs right
File.exists?(sample_file) ? sample_file = File.open(sample_file, 'r') : abort("Can't open the sample pool file!")

##### Class that stores information about each record from the sample key file
class Barcode_16s_record
  attr_accessor :pool, :barcode_num, :site_id, :patient, :sample
end

# Class that stores data about the reads
class Filtered_steps
attr_accessor :og_count, :lt_500bp, :gt_2000bp, :mapped, :singletons, :more_than_2_primers, :double_primers_and_correct, :not_primer_matched_by_usearch, :oriented, :singletons_retrieved, :percentage_retrieved			
end

#Argument (string) is the argument to check, arg_class is the class it should be, method is the method is being called from
def check_argument_type(argument, arg_class, method)
  abort("Argument #{argument} should be #{arg_class} but is not, in method #{method}") unless argument.class == arg_class
end



##### Method to write reads in fastq format
#fh = File handle, header = read header (string), sequence = read sequence, quality = phred quality scores
def write_to_fastq (fh, header, sequence, quality)

  fh       = file handle
  header   = string
  sequence = string
  quality  = array

	fh.write('@' + header + "\n")
	fh.write(sequence)
	fh.write("\n+\n")
	fh.write(quality + "\n")
end

##### Checking the number of projects based on the pacbio job id
def num_of_projects (sample_file)
	pb_projects = {}
	File.foreach(sample_file) do |entry|
  		next if $. == 1
  		arr = entry.split
  		pb_projects[arr.last] = [] unless pb_projects.has_key?(arr.last)
  		rec = Barcode_16s_record.new
  		rec.pool = arr[1]
  		rec.barcode_num = arr[2]
  		rec.site_id = arr[4]
  		rec.patient = arr[5]
  		rec.sample = arr[6]
  		pb_projects[arr.last].push(rec)
	end
	puts "Sanity check, number of projects are: " + pb_projects.count.to_s
	return pb_projects
end

##### Get the tarred_folder containging the fastq files and untar it 
def get_tarred_folder (id)
  	curr_file = "/data/pacbio/smrtanalysis_userdata/jobs/016/0" + id + "/data/barcoded-fastqs.tgz" if /^16/.match(id)
  	curr_file = "/data/pacbio/smrtanalysis_userdata/jobs/017/0" + id + "/data/barcoded-fastqs.tgz" if /^17/.match(id)
 	abort("The file #{curr_file} does not exist!") if !File.exists?(curr_file)
	`tar xvf #{curr_file}`
end

##### Get the ccs counts
def get_ccs_counts (id)   
	ccs_hash = {}                                                                                                               
	curr_files_ccs = "/data/pacbio/smrtanalysis_userdata/jobs/016/0" + id + "/data/*.ccs.h5" if /^16/.match(id)
  	curr_files_ccs = "/data/pacbio/smrtanalysis_userdata/jobs/017/0" + id + "/data/*.ccs.h5" if /^17/.match(id)
  	`python /data/shared/homes/archana/projects/ccs_passes.py #{curr_files_ccs} >> passes`
  	# Get a hash with the read header as key and ccs pass as value                                                                                             
  	passes_file = File.open("passes", "r")
  	passes_file.each do |line|
    		line = line.chomp
    		line_split = line.split("\s")
    		ccs_hash[line_split[0]] = line_split[1]
 	end
	return ccs_hash
end

##### Mapping reads to the human genome
def map_to_human_genome (base_name)                                                                                                               
      	#align all reads to the human genome                                                                                                                   
      	`bwa mem -t 15 ../human_g1k_v37.fasta corrected.fq > corrected.sam`

      	#sambamba converts sam to bam format                                                                                                                   
      	`sambamba view -S -f bam corrected.sam -o corrected.bam`

      	#Sort the bam file                                                                                                                                     
      	`sambamba sort -t15 -o corrected_sort.bam corrected.bam`

      	#filter the bam for only ‘not unmapped’ reads -> reads that are mapped                                                                                 
      	`sambamba view -F 'not unmapped' corrected.bam > #{base_name}_mapped.txt`
      	mapped_count = `cut -d ';' -f1 #{base_name}_mapped.txt| sort | uniq | wc -l`

      	#filter reads out for ‘unmapped’ -> we would use these for 16s id pipeline                                                                             
      	`sambamba view -F 'unmapped' corrected.bam > #{base_name}_unmapped.txt`

      	#convert the sam file to fastq                                                                                                                         
      	`grep -v ^@ #{base_name}_unmapped.txt | awk '{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}' > corrected_final.fastq`

	return mapped_count
end

##### Method for primer matching 
def primer_match (base_name)
	count_2_and_correct = 0                                                                                  
    count_1 = 0
    count_more_than_2 = 0
	record_hash = {}

	# Run usearch for primer matching
	if File.zero?("corrected_final.fastq")
		count_1, count_more_than_2, count_2_and_correct = 0
	else
      	`usearch -search_oligodb corrected_final.fastq -db /data/shared/homes/archana/projects/primers.fasta -strand both -userout #{base_name}_primer_map.txt -userfields query+target+qstrand+diffs+tlo+thi+qlo+qhi -matched #{base_name}_matched.fasta -notmatched #{base_name}_notmatched.fasta`

		# Opening the output from the primer matching step
    	primer_file = File.open("#{base_name}_primer_map.txt", "r")
      
    	# Storing the record header as the key and primer strand and sequence strand as the value. 
		primer_file.each do |line|
    		line_array = line.split("\t")
        	if record_hash.has_key?(line_array[0])
          		record_hash[line_array[0]].push([line_array[1], line_array[2]])
        	else
         		record_hash[line_array[0]] = [[line_array[1], line_array[2]]]
        	end
      	end
      	#puts record_hash.length                                                                                                                               

		record_hash.each do |key, value|
      		if value.size == 2
        		if value[0][0].eql?("reverse") && value[0][1].eql?("-") && value[1][0].eql?("forward") && value[1][1].eql?("+")
        			count_2_and_correct += 1
          		elsif value[0][0].eql?("reverse") && value[0][1].eql?("+") && value[1][0].eql?("forward") && value[1][1].eql?("-")
            		count_2_and_correct += 1
          		elsif value[0][0].eql?("forward") && value[0][1].eql?("+") && value[1][0].eql?("reverse") && value[1][1].eql?("-")
            		count_2_and_correct += 1
          		elsif value[0][0].eql?("forward") && value[0][1].eql?("-") && value[1][0].eql?("reverse") && value[1][1].eql?("+")
            		count_2_and_correct += 1
				end
   	  		elsif value.size == 1
          			count_1 += 1
        	elsif value.size > 2 
          			count_more_than_2 += 1
			end
		end 
	end 
    #puts to_complement_hash  
	return count_1, count_more_than_2, count_2_and_correct                                                                                                          
end	

def create_half_primer_files (primer_file_path)
	# Open the primer file
	primer_file = Bio::FlatFile.auto(primer_file_path, "r")

	# Get the database files with half forward and reverse from the primerfastafile
	out_fow_half = File.open("primer_half_fow.fasta", "w")
	out_rev_half = File.open("primer_half_rev.fasta", "w")

	primer_file.each do |entry|
		if entry.definition.include?("forward")
			half_len = entry.naseq.length/2
			# print to forward out
			out_fow_half.puts(">"+entry.definition)
			out_fow_half.puts(entry.naseq.upcase[half_len..-1])
			# print to reverse out
			out_rev_half.puts(">"+entry.definition)
			out_rev_half.puts(entry.naseq.upcase)
		elsif entry.definition.include?("reverse")
			rev_comple = entry.naseq.complement.upcase
			half_len = rev_comple.length/2
			# print to forward out
			out_fow_half.puts(">"+entry.definition)
			out_fow_half.puts(rev_comple)
			# print to reverse out
			out_rev_half.puts(">"+entry.definition)
			out_rev_half.puts(rev_comple[0..half_len])
		end
	end
	
	out_fow_half.close
	out_rev_half.close
	
	# Check to see if the files were created properly
	abort("!!!!The file primer_half_rev (required for primer matching with half the reverse primer sequence) does not exist!!!!") if !File.exists?("primer_half_rev.fasta")
	abort("!!!!The file primer_half_rev (required for primer matching with half the reverse primer sequence) is empty!!!!") if File.zero?("primer_half_rev.fasta")
	abort("!!!!The file primer_half_fow (required for primer matching with half the forward primer sequence) does not exist!!!!") if !File.exists?("primer_half_fow.fasta")
	abort("!!!!The file primer_half_fow (required for primer matching with half the forward primer sequence) is empty!!!!") if File.zero?("primer_half_fow.fasta")

end	

#### Process each file from the tarred folder 	
def process_each_file (samps, log, all_bc_reads, table, table_2, ccs_hash)
	samps.each do |rec|
    	log.puts("Pool: #{rec.pool} Barcode: #{rec.barcode_num}")
    	base_name = "#{rec.pool}_#{rec.barcode_num}_#{rec.site_id}_#{rec.patient}"
    	bc = "barcodelabel=#{base_name}\;"
    	fq = ''
    	current_bc = rec.barcode_num.to_i

		# Getting the appropriate fq file based on the barcode
    	if current_bc == 1
      		fq = '0001_Forward--0002_Forward.fastq'
      	elsif current_bc == 2
      		fq = "0003_Forward--0004_Forward.fastq"
      	elsif current_bc == 3
      		fq = "0005_Forward--0006_Forward.fastq"
      	elsif current_bc == 4
      		fq = "0007_Forward--0008_Forward.fastq"
      	end

		# Do the rest only if that fq file exits
    	if File.exists?(fq)
      		puts base_name
       	           
      		# Initialize log variable
      		filt = Filtered_steps.new
      		filt.og_count=0
      		filt.lt_500bp=0
      		filt.gt_2000bp=0
			filt.mapped=0
			filt.singletons=0
			filt.more_than_2_primers=0
			filt.double_primers_and_correct=0
			filt.not_primer_matched_by_usearch=0
            filt.oriented=0

      		# Add in the barcode to the fastq header
      		corrected = File.open('corrected.fq', 'w')
      		fh = Bio::FlatFile.auto(fq)

      		fh.each do |entry|
        		new_header = entry.definition.split[0]
        		if ccs_hash.has_key?(entry.definition.split[0])
          			ccs = "ccs=#{ccs_hash[entry.definition.split[0]]};"
          			new_header = new_header + ";" + bc + ccs
        		else
          			ccs = "ccs=nil;"
          			new_header = new_header + ";" + bc + ccs
        		end

        		write_to_fastq(corrected, "#{new_header}", entry.naseq.upcase, entry.quality_string) unless entry.naseq.size > 2000 || entry.naseq.size < 500
        		write_to_fastq(all_bc_reads, "#{new_header}", entry.naseq.upcase, entry.quality_string) unless entry.naseq.size > 2000 || entry.naseq.size < 500
        		filt.og_count  += 1
        		filt.lt_500bp  += 1 if entry.naseq.size < 500
        		filt.gt_2000bp += 1 if entry.naseq.size > 2000
			end

      		corrected.puts ""
      		corrected.close
      		File.delete(fq)

			# Calling the method which does human genome mapping
			mapped_count = map_to_human_genome (base_name) 
			filt.mapped = mapped_count.to_i

			# Calling the method for primer matching and getting statistics from it  
			count_1, count_more_than_2, count_2_and_correct = primer_match (base_name)
			filt.singletons = count_1
			filt.more_than_2_primers = count_more_than_2
			filt.double_primers_and_correct = count_2_and_correct

			# Dealing with the file which has sequences which are not primer matched by usearch                                                              
      		if File.exist?("#{base_name}_notmatched.fasta")
       			if File.zero?("#{base_name}_notmatched.fasta")
          			filt.not_primer_matched_by_usearch = 0
        		else
          			count_missing = `prinseq-lite -stats_info -fasta #{base_name}_notmatched.fasta`
          			count_missing_2 = count_missing.split("\t")[4].chomp.to_i
          			filt.not_primer_matched_by_usearch = count_missing_2
        		end
      		else
        		filt.not_primer_matched_by_usearch = 0
      		end	
				
			# Orient and get the number of oriented reads
			if File.exist?("corrected_final.fastq")
				if File.zero?("corrected_final.fastq")
					filt.oriented = 0
				else
      				`ruby parse_primer_matching_output_final.rb -f corrected_final.fastq -p #{base_name}_primer_map.txt -o #{base_name}_oriented_trimmed.fastq -s #{base_name}_singletons.fastq -m #{base_name}_singletons_primer_map.txt`
					if File.exist?("#{base_name}_oriented_trimmed.fastq")
						if File.zero?("#{base_name}_oriented_trimmed.fastq")
							filt.oriented = 0  
						else				
							count_oriented = `prinseq-lite -stats_info -fastq #{base_name}_oriented_trimmed.fastq`
      						count_oriented_2 = count_oriented.split("\t")[4].chomp.to_i
      						filt.oriented = count_oriented_2
						end
					end
				end
			else
				filt.oriented = 0
			end
			
			# Call the method which creates the primer files with half of the sequences
			create_half_primer_files("/data/shared/homes/archana/projects/primers.fasta")

			# Deals with singletons
			if File.exist?("#{base_name}_singletons.fastq")
				if File.zero?("#{base_name}_singletons.fastq")
					filt.singletons_retrieved = 0
					filt.percentage_retrieved = 0
				else
					count_forward_retrieved_2 = 0
					count_reverse_retrieved_2 = 0
					`ruby retrieve_singletons.rb -f #{base_name}_singletons.fastq -p #{base_name}_singletons_primer_map.txt -w #{base_name}_singletons_forward.fastq -r #{base_name}_singletons_reverse.fastq -x #{base_name}_singletons_primer_map_forward_trimmed.fastq -s #{base_name}_singletons_primer_map_reverse_trimmed.fastq`
					if File.exist?("#{base_name}_singletons_primer_map_forward_trimmed.fastq")			
						if File.zero?("#{base_name}_singletons_primer_map_forward_trimmed.fastq")
							count_forward_retrieved_2 = 0
						else
							count_forward_retrieved = `prinseq-lite -stats_info -fastq #{base_name}_singletons_primer_map_forward_trimmed.fastq`
							count_forward_retrieved_2 = count_forward_retrieved.split("\t")[4].chomp.to_i
						end
					end
					if File.exist?("#{base_name}_singletons_primer_map_reverse_trimmed.fastq")	
						if File.zero?("#{base_name}_singletons_primer_map_reverse_trimmed.fastq")
							count_reverse_retrieved_2 = 0
						else
							count_reverse_retrieved = `prinseq-lite -stats_info -fastq #{base_name}_singletons_primer_map_reverse_trimmed.fastq`
							count_reverse_retrieved_2 = count_reverse_retrieved.split("\t")[4].chomp.to_i
						end
					end
	      			filt.singletons_retrieved = count_forward_retrieved_2.to_i+count_reverse_retrieved_2.to_i
					filt.percentage_retrieved = (filt.singletons_retrieved.to_f/filt.singletons.to_f)*100
				end
			else	
				filt.singletons_retrieved = 0
				filt.percentage_retrieved = 0
			end	

			#puts "CHECK:", filt.percentage_retrieved
			#puts "SEQS FINALLY:", filt.oriented+filt.singletons_retrieved
	
      		# Write the filtered stats out to the log
      		table.puts("#{base_name}\t#{filt.og_count}\t#{filt.lt_500bp}\t#{filt.gt_2000bp}\t#{filt.mapped}\t#{filt.singletons}\t#{filt.more_than_2_primers}\t#{filt.double_primers_and_correct}\t#{filt.not_primer_matched_by_usearch}\t#{filt.oriented}\t#{filt.singletons_retrieved}\t#{filt.percentage_retrieved}")
			size_filt_total = filt.lt_500bp+filt.gt_2000bp	
			remains_after_size_filt = filt.og_count - size_filt_total
			remains_after_human_mapping = remains_after_size_filt - filt.mapped
			remains_after_primer_match_and_orienting= filt.oriented + filt.singletons_retrieved
			table_2.puts("#{base_name}\t#{filt.og_count}\t#{remains_after_size_filt}\t#{remains_after_human_mapping}\t#{remains_after_primer_match_and_orienting}")
		else
      		#log.puts("No file for sample #{id} barcode #{rec.barcode_num}")
      		table.puts("#{base_name}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0")
			table_2.puts("#{base_name}\t0\t0\t0\t0")
		end
	end 
end

##### Method which runs usearch and processes its output
def run_usearch_and_other (ee, uchime_db_file, utax_db_file, lineage_fasta_file) 
	# Merging all the files with oriented sequences
	`cat *_trimmed.fastq > all_ee#{ee}_oriented_primer_trimmed_merged.fastq`

	# Running the uparse pipeline
	`sh uparse_commands.sh all_ee#{ee} #{ee} #{uchime_db_file} #{utax_db_file}` 

	# Running the command to give a report of counts
	`ruby get_report.rb all_ee#{ee}`

	# Running blast on the OTUs                                                                                                                            
	#`usearch -usearch_global all_ee#{ee}_OTU_final.fasta -db ../16sMicrobial_ncbi_lineage.fasta -top_hit_only -id 0.9 -blast6out all_ee#{ee}_blast.txt -strand plus`
	`usearch -ublast all_ee#{ee}_OTU.fasta -db #{lineage_fasta_file} -top_hit_only -id 0.9 -blast6out all_ee#{ee}_blast.txt -strand both -evalue 0.01 -threads 15 -accel 0.3`

	# Running the script which maps between the blast file and the merged OTUs file                                                                         
	`ruby map_blast_with_otus.rb -u all_ee#{ee}_table_utax_map.txt -b all_ee#{ee}_blast.txt -o all_ee#{ee}_otu_blast_formatted.txt`

	# Running the script which dismantles the taxonomic assignments field
	`ruby parse_final_header.rb -h all_ee#{ee}_otu_blast_formatted.txt -o all_ee#{ee}_otu_blast_formatted2.txt`  
end

##### Main program calling all the methods

# Creating the files which are output by this script
# Lets make a log file for some stuff
log = File.open('log.txt', 'w')

# Let's also make a log for the counts
table = File.open("all_ee#{ee}_filter_sample_counts.txt", 'w')
table.puts("sample\tog_count\tlt_500bp\tgt_2000bp\tmapped\tsingletons\tmore_than_2_primers\tdouble_primers_which_passed\tnot_primer_matched_by_usearch\toriented\tsingletons_retrieved\tpercentage_singletons_retrieved")

table_2 = File.open("all_ee#{ee}_filter_sample_counts_summary.txt", 'w')
table_2.puts("sample\tog_count\tremains_after_size_filt\tremains_after_human_mapping\tremains_after_primer_matching_and_orienting")

# File with all original reads (>500 to <2000 in length)
all_bc_reads = File.open("all_ee#{ee}_bc_reads_size_filt.fq", 'w')

# Calling the method which counts the number of projects in the sample file
pb_projects = num_of_projects(sample_file)
#puts pb_projects

# Calling the methods which  
pb_projects.each do |id, samps|
	get_tarred_folder (id)
	ccs_hash = get_ccs_counts (id)
	#puts ccs_hash
	process_each_file(samps, log, all_bc_reads, table, table_2, ccs_hash)
end

# Calling the method which runs usearch and processes its output
run_usearch_and_other(ee, uchime_db_file, utax_db_file, lineage_fasta_file)

# Deleting unwanted files
File.delete('corrected.fq')
File.delete('corrected.bam')
File.delete('corrected.sam')
File.delete('corrected_sort.bam')
File.delete('corrected_final.fastq')
