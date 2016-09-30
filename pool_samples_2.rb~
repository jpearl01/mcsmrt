#!/usr/bin/env ruby
# -*- coding: utf-8 -*-

require 'bio'
require 'trollop'

#USAGE: ruby pool_samples.rb -s sample_key_BEI.txt -e 1 -c ../rdp_gold.fa -t ../16s_ncbi_database/16s_lineage_short_species_name_reference_database.udb -l ../16s_ncbi_database/16sMicrobial_ncbi_lineage.fasta

##### Input 
opts = Trollop::options do
  opt :samplefile, "File with all the sample information", :type => :string, :short => "-s"
  opt :eevalue, "Expected error at which you want filtering to take place", :type => :string, :short => "-e"
  opt :uchimedbfile, "Path to database file for the uchime command", :type => :string, :short => "-c"
  opt :utaxdbfile, "Path to database file for the utax command", :type => :string, :short => "-t"
  opt :lineagefastafile, "Path to FASTA file with lineage info for the ublast command", :type => :string, :short => "-l"
  opt :host_db, "Path to fasta file of host genome", :type => :string, :short => "-g"
end 

##### Assigning variables to the input and making sure we got all the inputs
opts[:samplefile].nil?       ==false  ? sample_file = opts[:samplefile]               : abort("Must supply a 'sample file': tab delimited file of sample information with '-s'")
opts[:eevalue].nil?          ==false  ? ee = opts[:eevalue]                           : abort("Must supply an Expected Error value with '-e'")
opts[:uchimedbfile].nil?     ==false  ? uchime_db_file = opts[:uchimedbfile]          : abort("Must supply a 'uchime database file' e.g. rdpgold.udb '-c'")
opts[:utaxdbfile].nil?       ==false  ? utax_db_file = opts[:utaxdbfile]              : abort("Must supply a 'utax database file' e.g. 16s_ncbi.udb '-t'")
opts[:lineagefastafile].nil? ==false  ? lineage_fasta_file = opts[:lineagefastafile]  : abort("Must supply a 'lineage fasta file' e.g. ncbi_lineage.fasta (for blast) '-l' ")
opts[:host_db].nil?          ==false  ? human_db = opts[:host_db]                     : abort("Must supply a fasta the host genome e.g. human_g1k.fasta  '-g' ")


##### Making sure we can open the sample file
File.exists?(sample_file) ? sample_file = File.open(sample_file, 'r') : abort("Can't open the sample pool file!")

##### Get the path to the directory in which the scripts exist 
script_directory = File.dirname(__FILE__)

##### Class that stores information about each record from the sample key file
class Barcode_16s_record
  attr_accessor :pool, :barcode_num, :site_id, :patient, :sample
end

#### Method to make sure the arguments given to a method are of the right type
#### Argument (string) is the argument to check, arg_class is the class it should be, method is the method is being called from
def check_argument_type(argument, arg_class, method)
  abort("Argument #{argument} should be #{arg_class} but is not, in method #{method}") unless argument.class == arg_class
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
	#puts pb_projects.inspect
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

#### Process each file from the tarred folder and get one file with all the reads
def process_each_file (samps, log, all_bc_reads, ccs_hash)
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
        
        write_to_fastq(all_bc_reads, "#{new_header}", entry.naseq.upcase, entry.quality_string) 
			end
		end
	end
end

# Work with the file which has all the reads
def process_all_bc_reads_file (all_bc_reads_to_process, all_reads_hash)
	all_bc_reads_to_process.each do |entry|
		# Fill the all_bc_hash with some basic info that we can get (read_name, barcode, ccs, length_pretrim)
		def_split = entry.definition.split(";")
		read_name = def_split[0]
		barcode = def_split[1].split("=")[1]
		ccs = def_split[2].split("=")[1].to_i
		all_reads_hash[read_name] = Array.new(12)
		all_reads_hash[read_name][0] = read_name
		all_reads_hash[read_name][1] = ccs
		all_reads_hash[read_name][2] = barcode
		all_reads_hash[read_name][6] = entry.naseq.size

		# Get ee_pretrim
		# Call the get_ee_from_fq_file method

	end	
end

# Make this a generic method whcih takes an fq file as argument and returns a hash with the read name and ee
def get_ee_from_fq_file (file)
	`usearch -fastq_filter #{file} -fastqout #{} -fastq_maxee #{ee} -fastq_qmax 75 -fastq_eeout -sample all`
end

# Creating the files which are output by this script
# Lets make a log file for some stuff
log = File.open('log.txt', 'w')

# File with all original reads
all_bc_reads = File.open("all_ee#{ee}_bc.fq", 'w')

# Calling the method which counts the number of projects in the sample file
pb_projects = num_of_projects(sample_file)
#puts pb_projects

# Calling the methods which help in producing one file with all the reads
pb_projects.each do |id, samps|
  get_tarred_folder (id)
  ccs_hash = get_ccs_counts (id)
  #puts ccs_hash
  process_each_file(samps, log, all_bc_reads, ccs_hash)
end

all_reads_hash = {}
all_bc_reads_to_process = Bio::FlatFile.auto("all_ee#{ee}_bc.fq")
process_all_bc_reads_file(all_bc_reads_to_process, all_reads_hash)
puts all_reads_hash




