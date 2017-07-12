# MC-SMRT
## Microbiome Classifier for SMRT PacBio data

### Introduction:
MCSMRT is a tool for microbiome analysis of PacBio data. This tool is a pipeline to go from raw PacBio data to clustered sequences (OTUs or Operational Taxonomic Units).  The outputs include a table of the number of reads assigned to each OTU and taxonomic lineage associated with each OTU. Another output consists of general information regarding each read, for example, CCS count, expected error, length, primer matching result, etc. 

### Installation and Dependencies: 
1. Ruby v2.2.1 or greater
   To install dependencies run the command 
   `$ bundle`
2. BWA https://sourceforge.net/projects/bio-bwa/files/
3. Sambamba http://lomereiter.github.io/sambamba/
4. Usearch v8.1 from http://www.drive5.com/usearch/download.html .  
   `$ln -s usearch8.1.*_i86linux* ~/bin/usearch`

### Requirements:
Fasta file of forward and reverse primer sequences. (For more information on FASTA file formats, refer https://en.wikipedia.org/wiki/FASTA_format).

Fasta formatted taxonomy classification database (NCBI database included in this repository)

Usearch formatted taxonomy classification database see: If you want to create your own udb file, refer http://www.drive5.com/usearch/manual/cmd_makeudb_utax.html (NCBI databse included repository)

rdp_gold database for chimera detection: http://drive5.com/uchime/gold.fa

Host genome file, for filtering sequences mapping to host genome

### Usage: 
First and foremost, you will need a FASTA file with primer sequences used for your project. Go ahead and create the file with all the primer sequences used and store it in your working directory.
This file is used for matching sequences to primers and trimming (which is optional). This file is referred to as PRIMERS_DB and is one of the inputs for the mcsmrt.rb script.

There are 2 methods/routes you can choose from, in order to run your PacBio data through the microbiome classifier. 

#### Method 1: For cases when you have already demultiplexed your samples
  * Step 1: Run get_fastqs.rb in order to obtain CCS counts and barcode information as a part of the FASTQ headers. This step also copies the demultiplexed FASTQ files into a folder. Folder name in which the files are stored should be given by you as an argument.  
    Command:    
    ruby get_fastqs.rb [-h] [-s SAMPLE_INFO_FILE] [-o OUTPUT_FOLDER_NAME]  

  * Step 2: Run mcsmrt.rb in order to obtain detailed information about each read and find the clusters (OTUs) created. This step in-turn can be used in 2 ways. 
    Specify that all the files in a given directory can be used for clustering (usually when all the files belong to the same project), or, provide a file with a list of file names which can be clustered together.
    Command:  
    ruby mcsmrt.rb [-h]  
[-a] / [-i LIST_OF_FILES_FOR_CLUSTEIRNG] [-f FOLDER_NAME] [-m TRIMMING]  
[-e EXPECTED_ERROR] [-s CCS_COUNT] [-x MAXIMUM_LENGTH] [-n MINIMUM_LENGTH]  
[-c UCHIME_DB] [-t UTAX_DB] [-l BLAST_DB] [-g HOST_GENOME_DB] [-p PRIMERS_DB]                                                                                        

#### Method 2: For cases when you have to demultiplex your samples and then run it through the microbiome classifier
  * Step 1: Demultiplex using Rachel’s demultiplexing pipeline (https://github.com/rehrlich/ccs_smrt_pipe). This pipeline also results in FASTQ files which have a format that is compatible with mcsmrt.rb script for microbiome analysis. 
  * Step 2: Run mcsmrt.rb as described above. 

### Arguments explained:

#### Arguments for get_fastqs.rb:
  * SAMPLE_INFO_FILE (-s) – This is the file which will have a list of all the PacBio jobs which are demutiplexed. The FASTQ files of the jobs given in this list will be extracted and the headers of each FASTQ sequence will be added with additional tags like CCS count and sample information. The header of this file (first row) should have column names corresponding to PB_jobid, data_path, forward_barcode, reverse_barcode and sample_name. These column names HAVE TO BE exactly as is described here because the program initializes data in each column based on these column names. Data in each column is described as follows:  
      1. PB_jobid – The ID number assigned to each job when the data is demultiplexed on SMRT portal.    
      2. data_path – Path to where the demultiplexed data for each job is located.   
      3. forward_barcode – Name of the forward primer used for the samples in the respective job, as given in the PRIMERS_DB file.  
      4. reverse_barcode - Name of the reverse primer used for the samples in the respective job, as given in the PRIMERS_DB file.  
      5. sample_name – This is the name given to each sample. This is the one that is going to be added in the FASTQ sequence header with a tag of “barcodelabel”. So, if you want any information to be kept track of, add it as a sample name. Multiple things can be kept track of in the sample name, separated by a “_”. For example, if I want to keep track of patient ID and sample ID in this location, give it the sample name “Pat123_Samp167” where Pat123 corresponds to the patient ID and Samp167 corresponds to the sample ID. This way all this information will be associated with each sequence and can later be tracked easily. Mandatorily, each forward and reverse barcode pair is given a unique number and this number is also added in the beginning of the “barcodelabel” tag in the fastq file.    
  * OUTPUT_FOLDER_NAME (-o) – Name of the folder in which the FASTQ files (with ccs and sample information in the header) are going to be stored. The name of each FASTQ file in this folder is going to be the same as the sample name given in the SAMPLE_INFO_FILE.   

#### Arguments for mcsmrt.rb:
  * ALL (-a) / LIST_OF_FILES_FOR_CLUSTEIRNG (-i) – If you want to specify particular files for clustering, you should create another file with a list of all the FASTQs you want to concatenate for clustering. This file must have one file name in each row.  If you want all the files in the directory to be clustered together, use the –a option. 
  * FOLDER_NAME (-f) – Name of the folder (preferably the full path), in which the files for microbiome analysis and/or clustering exist. So, the files you provided for this script with the -a/-i option should all be stored in the directory given with this argument. 
  * EXPECTED_ERROR (-e), Default (1.0) – Maximum expected error above which sequences are filtered out (refer  http://www.drive5.com/usearch/manual/expected_errors.html for more information on expected error).
  * CCS_COUNT (-s), Default (5) – Minimum CCS count below which sequences are filtered out. 
  * MAXIMUM_LENGTH (-x), Default (2000) – Maximum length above which sequences are filtered out. 
  * MINIMUM_LENGTH (-n), Default (500) – Minimum length below which sequences are filtered out. 
  * UCHIME_DB (-c) – The database file which has high quality 16s sequences. I use the rdp_gold database for this purpose. This file is used to filter out OTU sequences which are chimeras. 
  * UTAX_DB (-t) – The udb format database file which has a lineage assigned to each 16s sequence. These sequences are also trained for getting the confidence at each taxonomic level. This is the file that is used to assign taxonomy and confidences to OTU sequences.
We used our custom made database file which uses a tool called Lineanator for this purpose (https://github.com/bhatarchanas/lineanator). Lineanator uses NCBI 16s database and NCBI taxonomy database to create a udb format file for assigning taxonomy. If you want to create your own udb file, refer http://www.drive5.com/usearch/manual/cmd_makeudb_utax.html. 
  * BLAST_DB (-l) – This is the FASTA file that was used to create the udb format file. In this script, the file is used to get a strain name, alignment length and percent identity for each OTU sequence. 
  * HOST_GENOME_DB (-g) – This file is used as a reference to filter out those sequences which map to the host genome. It should be a FASTA format file. 
  * PRIMERS_DB (-p) – This is the FASTA format primer file that you created for primer matching and trimming.
  * TRIMMING (-m), Default (yes) – Give “yes” if you want the sequences to be trimmed/primer sequences to be removed, and no if you don’t. 

### Built with:  
ruby 2.2.1p85 (2015-02-26 revision 49769) [x86_64-linux]
