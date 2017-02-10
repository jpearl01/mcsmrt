require 'bio'
require 'trollop'

#USAGE: ruby get_fastqs.rb -s #{sample_key_file} -o #{all_bc_reads_output_folder} 

##### Input 
opts = Trollop::options do
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

##### Checking the number of projects based on the pacbio job id
def num_of_projects (sample_file)
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
    		rec = Barcode_16s_record.new
    		rec.fow_barcode = entry_split[fow_barcode_ind]
    		rec.rev_barcode = entry_split[rev_barcode_ind]
    		rec.sample = entry_split[sample_ind]
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
	if samps[0].end_with?("/")
		curr_file = samps[0] + "data/barcoded-fastqs.tgz"
		abort("The file #{curr_file} does not exist!") if !File.exists?(curr_file)
		`tar xvf #{curr_file} -C raw_data/`
	else
		curr_file = samps[0] + "/data/barcoded-fastqs.tgz"
		abort("The file #{curr_file} does not exist!") if !File.exists?(curr_file)
		`tar xvf #{curr_file} -C raw_data/`
	end
end

##### Get the ccs counts
def get_ccs_counts (samps, script_directory)   
	if samps[0].end_with?("/")                                                                                                      
  		curr_files_ccs = samps[0] + "data/*.ccs.h5"
  		`python #{script_directory}/ccs_passes.py #{curr_files_ccs} >> passes`
  	else
		curr_files_ccs = samps[0] + "/data/*.ccs.h5"
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

  	File.delete(passes_file)
  	return ccs_hash
end

##### Get barcode nums or nick names for each barcode pair 
def get_barcode_nums (samps)
	barcode_hash = {}
	count = 0
	samps[1].each do |rec|
		#puts rec.rev_barcode
		#puts rec.fow_barcode
		key = "#{rec.fow_barcode}--#{rec.rev_barcode}"
		if barcode_hash.key?(key)
			next
		else
			count += 1
			barcode_hash[key] = count
		end 
	end

	barcode_num_out_file = File.open("barcode_num_for_pairs.txt", "w")
	barcode_num_out_file.puts("barcode_pair\tbarcode_num")
	barcode_hash.each do |key, value|
		barcode_num_out_file.puts(key+"\t"+value.to_s)
	end

	return barcode_hash
end

##### Process each file from the tarred folder and get one file with all the reads
def process_each_file (samps, log, output_folder, ccs_hash)
  	# Write the status of whether the file was processed or not in the log file
  	log.puts("STATUS: BarcodeNum_Sample")

  	# Call the method which create a barcode nick-names file
  	barcode_hash = get_barcode_nums(samps)

  	# Loop through each sample pointing to a job id
  	samps[1].each do |rec|
  		# Opening the file with the original reads are written along with the extra info
  		all_bc_reads = File.open("#{output_folder}/#{rec.sample}.fq", "w")	

  		barcode_num = barcode_hash["#{rec.fow_barcode}--#{rec.rev_barcode}"]
    	base_name = "#{barcode_num}_#{rec.sample}"
    	bc = "barcodelabel=#{base_name}\;"
    	fq = ''
    
    	# Getting the appropriate fq file based on the barcode
      	fq = "raw_data/#{rec.fow_barcode}--#{rec.rev_barcode}.fastq"

    	# Do the rest only if that fq file exits
    	if File.exists?(fq)
    		log.puts("FOUND: #{base_name}")
      		puts base_name
      
      		# Add in the barcode and ccs count to the fastq header
      		fh = Bio::FlatFile.auto(fq)			
      		fh.each do |entry|
				new_header_raw = entry.definition.split[0]
        		if ccs_hash.has_key?(entry.definition.split[0])
          			ccs = "ccs=#{ccs_hash[entry.definition.split[0]]};"
          			new_header = new_header_raw + ";" + bc + ccs
        		else
          			ccs = "ccs=nil;"
          			new_header = new_header_raw + ";" + bc + ccs
        		end
        
       			# Write all the reads into the output file which was provided
        		write_to_fastq(all_bc_reads, "#{new_header}", entry.naseq.upcase, entry.quality_string) 
			end
		else
			log.puts("NOT FOUND: #{base_name}")	
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

# Calling the method which counts the number of projects in the sample file
pb_projects = num_of_projects(sample_file)
#puts pb_projects

# Calling the methods which help in producing files with all the reads in the output folder provided.
pb_projects.each do |id, samps|
	#puts id, samps
	
	# Create a directory which will hold the raw data fq files 
	if Dir.exists?("raw_data")
		Dir.glob("raw_data/*.fq") do |fq_file|
  			fq_file.delete
		end
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
`rm -rf raw_data`

##########################################################