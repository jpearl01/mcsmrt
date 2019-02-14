#!/usr/bin/env ruby

require 'bio'
require 'optimist'
require 'fileutils'

#USAGE: ruby ../mcsmrt/get_fastqs.rb -s #{sample_key_file} -o #{all_bc_reads_output_folder} 

##### Input 
opts = Optimist::options do
	opt :samplefile, "File with all the sample information", :type => :string, :short => "-s"
	opt :outfolder, "Output FASTQ file which has all the reads in it.", :type => :string, :short => "-o"
end 

##### Assigning variables to the input and making sure we got all the inputs
opts[:samplefile].nil? ==false  ? sample_file = opts[:samplefile]   : abort("Must supply a 'sample file' which is a tab delimited file of sample information with '-s'")
opts[:outfolder].nil?  ==false  ? output_folder = opts[:outfolder]  : abort("Must supply an 'output file name' which is the name of the fq file with all the reads from the jobs given in the sample sheet with '-o'")

##### Making sure we can open the sample file
File.exists?(sample_file) ? sample_file = File.open(sample_file, 'r') : abort("Can't open the sample pool file that is given to the '-s' argument!")

##### Get the path to the directory in which the scripts exist 
script_directory = File.dirname(__FILE__)

##### Class that stores information about each record from the sample key file
class Barcode_16s_record
  	attr_accessor :fow_barcode, :rev_barcode, :sample, :path
end

##### Method to write reads in fastq format
#fh = File handle, header = read header (string), sequence = read sequence (string), quality = phred quality scores (array)
def write_to_fastq (fh, header, sequence, quality)
  	fh.write('@' + header + "\n")
  	fh.write(sequence)
  	fh.write("\n+\n")
  	fh.write(quality + "\n")
end

def prompt_to_change_sample_names (sample_file)
	prompt_hash = {}
	sample_ind = ""
	special = "?!<>',?[]}{=-)(*&^%$#`~{}"
	regex = /[#{special.gsub(/./){|char| "\\#{char}"}}]/

	sample_file.each_with_index do |entry, ind|
		if ind == 0
			header_split = entry.chomp.split("\t")
			sample_ind = header_split.index("sample_name")
		else
			entry_split = entry.chomp.split("\t")
			sample = entry_split[sample_ind]
			if entry_split[sample_ind][0] =~ /^[0-9].*/ or entry_split[sample_ind] =~ regex
				prompt_hash[sample] = "yes"
			else
				prompt_hash[sample] = "no"
			end
		end
	end

	if prompt_hash.has_value?("yes")
		puts "Sample names need to be changed. Should I go ahead and change it myself? Answer in yes or no."
		answer = gets.chomp
		if answer == "yes"
			puts "Slightly modifying sample names and moving ahead...."
		elsif answer == "no"
    		abort("Correct the sample names manually and re-run this script. Use samples names without special characters, underscores are allowed. Sample names should always begin with an alphabet.")
    	else
    		abort("You did not answer in yes or no. Aborting the code now. Run the script again.")
   		end
   	end

   	sample_file.close()
	return prompt_hash

end

##### Checking the number of projects based on the pacbio job id
def num_of_projects (sample_file, prompt_hash)
  	pb_projects = {} # A hash which points to one array for each job_id. Structure is pb_projects[pb_jobid] = [path_to_files, [sample1, sample2]]. Each sample is an object of the Barcode_16s_record class!

  	fow_barcode_ind = 0
  	rev_barcode_ind = 0
	site_id_ind = 0
	patient_ind = 0
	sample_ind = 0
	pb_jobid_ind = 0
	data_path_ind = 0
	
	sample_file.each_with_index do |entry, ind|
		if ind == 0
			header_split = entry.chomp.split("\t")
			fow_barcode_ind = header_split.index("forward_barcode")
			rev_barcode_ind = header_split.index("reverse_barcode")
			sample_ind = header_split.index("sample_name")
			pb_jobid_ind = header_split.index("PB_jobid")
			data_path_ind = header_split.index("data_path")
			#puts pb_jobid_ind, pool_ind, barcode_num_ind
		else
			entry_split = entry.chomp.split("\t")
			pb_projects[entry_split[pb_jobid_ind]] = ["", []] unless pb_projects.has_key?(entry_split[pb_jobid_ind])

			special = "?!<>',?[]}{=-)(*&^%$#`~{}"
			regex = /[#{special.gsub(/./){|char| "\\#{char}"}}]/

			# Calling the method which checks sample names
			if prompt_hash[entry_split[sample_ind]] == "yes"
				sample_name_mod = ""
				if entry_split[sample_ind][0] =~ /^[0-9].*/
					sample_name_mod = "Samp_#{entry_split[sample_ind]}"
					#puts sample_name_mod
					if sample_name_mod =~ regex
						sample_name_mod = sample_name_mod.gsub(/[^\w]/, '_')
						#puts sample_name_mod
					end
				end

				# Populating the rec class and pb_projects hash
    			rec = Barcode_16s_record.new
    			rec.fow_barcode = entry_split[fow_barcode_ind]
    			rec.rev_barcode = entry_split[rev_barcode_ind]
    			rec.sample = sample_name_mod
    		else
    			rec = Barcode_16s_record.new
    			rec.fow_barcode = entry_split[fow_barcode_ind]
    			rec.rev_barcode = entry_split[rev_barcode_ind]
    			rec.sample = entry_split[sample_ind]
			end

    		pb_projects[entry_split[pb_jobid_ind]][0] = entry_split[data_path_ind]
    		pb_projects[entry_split[pb_jobid_ind]][1].push(rec)

		end
	end

  	puts "Sanity check, number of projects are: " + pb_projects.count.to_s
	#puts pb_projects.inspect
  	return pb_projects
end

##### Get the tarred_folder containging the fastq files and untar it 
def get_tarred_folder (samps)
	#puts samps.inspect
	if samps[0].end_with?("/")
		curr_file = samps[0] + "data/barcoded-fastqs.tgz"
		abort("The file #{curr_file} does not exist!") if !File.exists?(curr_file)
		`tar xvf #{curr_file} -C raw_data/`
	else
		curr_file = "#{samps[0]}/" + "data/barcoded-fastqs.tgz"
		abort("The file #{curr_file} does not exist!") if !File.exists?(curr_file)
		`tar xvf #{curr_file} -C raw_data/`
	end
end

##### Get the ccs counts
def get_ccs_counts (samps, script_directory) 

	if File.exists?("passes")
		FileUtils.rm_f("passes")
	end

	if samps[0].end_with?("/")                                                                                                      
  		curr_files_ccs = samps[0] + "data/*.ccs.h5"
  		`python #{script_directory}/ccs_passes.py #{curr_files_ccs} >> passes`
  	else
		curr_files_ccs = "#{samps[0]}/" + "data/*.ccs.h5"
		`python #{script_directory}/ccs_passes.py #{curr_files_ccs} >> passes`
  	end
	
	# Get a hash with the read header as key and ccs pass as value        
	ccs_hash = {}                                                                                               
  	passes_file = File.open("passes", "r")
  	passes_file.each do |line|
    	line = line.chomp
    	line_split = line.split("\s")
    	ccs_hash[line_split[0]] = line_split[1]
  	end

  	#File.delete(passes_file)
  	return ccs_hash
end

##### Process each file from the tarred folder and get one file with all the reads
def process_each_file (samps, log, output_folder, ccs_hash)
  	# Write the status of whether the file was processed or not in the log file
  	log.puts("STYLE: Sample;Barcode")

  	# Loop through each sample pointing to a job id
  	samps[1].each do |rec|
  		#puts samps[1].inspect

  		# Opening the file with the original reads are written along with the extra info
  		if File.exists?("#{output_folder}/#{rec.sample}.fq")
  			abort("Unique sample names must be provided. Sample name #{rec.sample} exists multiple times.")
  		else
  			all_bc_reads = File.open("#{output_folder}/#{rec.sample}.fq", "w")
  		end

    	bc = "barcodelabel=#{rec.sample}\;barcode=#{rec.fow_barcode}--#{rec.rev_barcode}\;"
    	fq = ''
    
    	# Getting the appropriate fq file based on the barcode
      	fq = "raw_data/#{rec.fow_barcode}--#{rec.rev_barcode}.fastq"

    	# Do the rest only if that fq file exits
    	if File.exists?(fq)
    		log.puts("FOUND: #{rec.sample};#{rec.fow_barcode}--#{rec.rev_barcode}")
      		puts "#{rec.sample};#{rec.fow_barcode}--#{rec.rev_barcode}"
      
      		# Add in the barcode and ccs count to the fastq header
      		fh = Bio::FlatFile.auto(fq)			
      		fh.each do |entry|
				new_header_raw = entry.definition.split[0]
        		if ccs_hash.has_key?(entry.definition.split[0])
          			ccs = "ccs=#{ccs_hash[entry.definition.split[0]]};"
          			new_header = new_header_raw + ";" + ccs + bc
        		else
          			ccs = "ccs=nil;"
          			new_header = new_header_raw + ";" + ccs + bc
        		end
        
       			# Write all the reads into the output file which was provided
        		write_to_fastq(all_bc_reads, "#{new_header}", entry.naseq.upcase, entry.quality_string) 
			end
		else
			log.puts("NOT FOUND: #{rec.sample};#{rec.fow_barcode}--#{rec.rev_barcode}")	
		end

		# Close the file in which all the reads were written
		all_bc_reads.close
	end
end

##################### MAIN PROGRAM #######################

# Lets make a log file for some stuff
log = File.open('log.txt', 'w')

# Creating a folder with the name of the given output directory
if Dir.exists?("#{output_folder}")
	abort("Output folder already exists. Either delete the folder if it is not needed or enter a new output folder name and re-run the command!")
else
	Dir.mkdir("#{output_folder}")
end

# Calling the method which checks sample names and makes sure the user is ok with modying them within the program
prompt_hash = prompt_to_change_sample_names(sample_file)

# re-opening the sample file since it was read and closed in the promt method
sample_file = File.open(sample_file, 'r')
# Calling the method which counts the number of projects in the sample file
pb_projects = num_of_projects(sample_file, prompt_hash)
puts pb_projects.inspect

# Calling the methods which help in producing files with all the reads in the output folder provided.
pb_projects.each do |id, samps|
	#puts id, samps
	
	# Create a directory which will hold the raw data fq files 
	if Dir.exists?("raw_data")
		FileUtils.rm_rf("raw_data")
		Dir.mkdir("raw_data")
	else
		Dir.mkdir("raw_data")
	end

	# Getting the raw data files and processing them
	get_tarred_folder(samps)
	ccs_hash = get_ccs_counts(samps, script_directory)
	#puts ccs_hash
  	
	# Calling the method which produces the final folder with the fastq files
  	process_each_file(samps, log, output_folder, ccs_hash)
end

# Delete the raw data files 
FileUtils.rm_rf("raw_data")

##########################################################
