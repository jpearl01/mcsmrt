[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [](#lang-us) ![ruby in bioinformatics ftw](https://img.shields.io/badge/Language-ruby-steelblue.svg)


# MC-SMRT
## Microbiome Classifier for SMRT PacBio data

### Introduction:
MCSMRT is a tool for microbiome analysis of PacBio FL16S sequence data. This pipeline is intended to help a user go from raw PacBio data to taxonomically classified OTU (Operational Taxonomic Units) sequence clusters. The outputs of greatest interest to users are likely a table of read counts assigned to each OTU and its corresponding taxonomic lineage, and a file consisting of specific information regarding each read, (e.g., CCS count, expected error, length, primer matching result, etc.). 

### Installation and Dependencies: 
1. Sample Data, and classification database can be [downloaded here](https://drive.google.com/open?id=1UJZBU3PhEVq8lUGcjPcs2s2LbqjsQctA)

2. Install Ruby v2.2.1 or greater 
   For Fedora:
   ```
   $ sudo dnf install ruby
   ```

   For Centos/RHEL
   ```
   $ sudo yum install ruby
   ```

   Install the 'bundler' gem:
   ```
   $ gem install bundler
   ```

   To install dependencies run the command   
   ```
   $ bundle
   ```

3. Install [BWA](https://sourceforge.net/projects/bio-bwa/files/)
   Follow the directions to compile the bwa exectuable, then add it to your path. The easiest way:
   ```
   $ ln -s /path/to/bwa ~/bin/
   ```

4. Install [Sambamba](http://lomereiter.github.io/sambamba/)
   Download precompiled binary, unzip and add to your path, for example:
   ```
   $ tar xvf sambamba* && mv sambamba ~/bin/
   ```


5. Install [Usearch v8.1](http://www.drive5.com/usearch/download.html). Download the 32bit version after agreeing to the license. After downloading make sure to make the file exectuable:
```
$ chmod +x usearch8.1.*
```

Add a softlink to your path (making sure the executable name is 'usearch'). For example:
   ```
   $ln -s /home/user/apps/usearch8.1.*_i86linux* ~/bin/usearch
   ```

### Data Prerequisites:
[Fasta file](https://en.wikipedia.org/wiki/FASTA_format) of forward and reverse primer sequences used in your specific PCR protocol.

USEARCH formatted taxonomy classification database is required for assigning taxonomy to OTUs. If you want to format your own udb file, refer [see here](http://www.drive5.com/usearch/manual/cmd_makeudb_utax.html) (our formatted NCBI FL16S databse is available [here](https://drive.google.com/open?id=1UJZBU3PhEVq8lUGcjPcs2s2LbqjsQctA)). This database was created using a curated set of the NCBI's 16s repository and taxonomy data. The tool which generates this database is called Lineanator and is available [here](https://github.com/bhatarchanas/lineanator). 

Fasta formatted [taxonomy classification database](https://drive.google.com/open?id=1UJZBU3PhEVq8lUGcjPcs2s2LbqjsQctA). Same as above, except formatted as a FASTA file. This file can also be created by a user with [Lineanator](https://github.com/bhatarchanas/lineanator).  

Clustering results from the database that was used for assigning taxonomy.

[RDP gold database](http://drive5.com/uchime/gold.fa) for chimera detection 

Host genome file, for filtering sequences mapping to host genome.  This may or may not make sense for you depending on the provenance of your samples. As this tool was created using primarily clinical samples, filtering out human DNA was often an important quality control step.

### Usage: 
First you will need a FASTA file with the PCR primer sequences used for your 16s amplification. Create or acquire the file with all the primer sequences used and store it in your working directory.
This file is used for matching sequences to primers and trimming (which is optional). This file is referred to as PRIMERS_DB and is one of the inputs for the mcsmrt.rb script.

There are 2 methods/routes you can use to run your PacBio data through this microbiome classifier. 

#### Method 1: For cases when you have already demultiplexed your samples
  * Step 1: Run get_fastqs.rb in order to insert CCS counts, and barcode information into the FASTQ read headers. This step also copies the demultiplexed FASTQ files into a user defined output folder. You will need both your fastq files, and a 'sample_info' file (format described below).
    Command:    
    `ruby get_fastqs.rb [-h] [-s SAMPLE_INFO_FILE] [-o OUTPUT_FOLDER_NAME]`

  * Step 2: Run mcsmrt_v1.rb in order to obtain detailed information about each read and find the clusters (OTUs) created. This step in-turn can be used in 2 ways. 
    Specify that all the files in a given directory can be used for clustering (usually when all the files belong to the same project), or, provide a file with a list of file names which can be clustered together.
    Command:  
  ```
    ruby mcsmrt_v1.rb [-h]  
	[-a] / [-i LIST_OF_FILES_FOR_CLUSTEIRNG] [-f FOLDER_NAME] [-m TRIMMING]  
	[-e EXPECTED_ERROR] [-s CCS_COUNT] [-x MAXIMUM_LENGTH] [-n MINIMUM_LENGTH]  
	[-c UCHIME_DB] [-t UTAX_DB] [-l BLAST_DB] [-g HOST_GENOME_DB] [-p PRIMERS_DB]
	[-d THREADS] [-b NCBI_CLUSTERED_FILE] [-v VERBOSE]                                                                                        
```

#### Method 2: For cases when you have to demultiplex your samples and then run it through the microbiome classifier
  * Step 1: Demultiplex using [this demultiplexing pipeline](https://github.com/rehrlich/ccs_smrt_pipe). This pipeline also results in FASTQ files which have a format that is compatible with mcsmrt.rb script for microbiome analysis. 
  * Step 2: Run mcsmrt.rb as described above. 

### Arguments explained:

#### Arguments for get_fastqs.rb:
  * `SAMPLE_INFO_FILE (-s)` – This is the file which will have a list of all the PacBio jobs which are demutiplexed. The FASTQ files of the jobs given in this list will be extracted and the headers of each FASTQ sequence will be added with additional tags like CCS count and sample information. The header of this file (first row) should have column names corresponding to PB_jobid, data_path, forward_barcode, reverse_barcode and sample_name. These column names HAVE TO BE exactly as is described here because the program initializes data in each column based on these column names. Data in each column is described as follows:  
      1. PB_jobid – The ID number assigned to each job when the data is demultiplexed on SMRT portal.    
      2. data_path – Path to where the demultiplexed data for each job is located.   
      3. forward_barcode – Name of the forward primer used for the samples in the respective job, as given in the PRIMERS_DB file.  
      4. reverse_barcode - Name of the reverse primer used for the samples in the respective job, as given in the PRIMERS_DB file.  
      5. sample_name – This is the name given to each sample. This is the one that is going to be added in the FASTQ sequence header with a tag of “barcodelabel”. So, if you want any information to be kept track of, add it as a sample name. Multiple things can be kept track of in the sample name, separated by a “_”. For example, if I want to keep track of patient ID and sample ID in this location, give it the sample name “Pat123_Samp167” where Pat123 corresponds to the patient ID and Samp167 corresponds to the sample ID. This way all this information will be associated with each sequence and can later be tracked easily. Mandatorily, each forward and reverse barcode pair is given a unique number and this number is also added in the beginning of the “barcodelabel” tag in the fastq file.    
  * OUTPUT_FOLDER_NAME (-o) – Name of the folder in which the FASTQ files (with ccs and sample information in the header) are going to be stored. The name of each FASTQ file in this folder is going to be the same as the sample name given in the SAMPLE_INFO_FILE.   

#### Arguments for mcsmrt.rb:
  * `ALL (-a) / LIST_OF_FILES_FOR_CLUSTEIRNG (-i)` – If you want to specify particular files for clustering, you should create another file with a list of all the FASTQs you want to concatenate for clustering. This file must have one file name in each row.  If you want all the files in the directory to be clustered together, use the –a option. 
  * `FOLDER_NAME (-f)` – Name of the folder (preferably the full path), in which the files for microbiome analysis and/or clustering exist. So, the files you provided for this script with the -a/-i option should all be stored in the directory given with this argument. 
  * `EXPECTED_ERROR (-e)`, Default (1.0) – Maximum expected error above which sequences are filtered out (refer  http://www.drive5.com/usearch/manual/expected_errors.html for more information on expected error).
  * `CCS_COUNT (-s)`, Default (5) – Minimum CCS count below which sequences are filtered out. 
  * `MAXIMUM_LENGTH (-x)`, Default (2000) – Maximum length above which sequences are filtered out. 
  * `MINIMUM_LENGTH (-n)`, Default (500) – Minimum length below which sequences are filtered out. 
  * `UCHIME_DB (-c)` – The database file which has high quality 16s sequences. I use the rdp_gold database for this purpose. This file is used to filter out OTU sequences which are chimeras. 
  * `UTAX_DB (-t)` – The udb format database file which has a lineage assigned to each 16s sequence. These sequences are also trained for getting the confidence at each taxonomic level. This is the file that is used to assign taxonomy and confidences to OTU sequences.
We used our custom made database file which uses a tool called Lineanator for this purpose (https://github.com/bhatarchanas/lineanator). Lineanator uses NCBI 16s database and NCBI taxonomy database to create a udb format file for assigning taxonomy. If you want to create your own udb file, refer http://www.drive5.com/usearch/manual/cmd_makeudb_utax.html. 
  * `BLAST_DB (-l)` – This is the FASTA file that was used to create the udb format file. In this script, the file is used to get a strain name, alignment length and percent identity for each OTU sequence. 
  * `HOST_GENOME_DB (-g)` – This file is used as a reference to filter out those sequences which map to the host genome. It should be a FASTA format file. 
  * `PRIMERS_DB (-p)` – This is the FASTA format primer file that you created for primer matching and trimming.
  * `TRIMMING (-m)`, Default (yes) – Give “yes” if you want the sequences to be trimmed/primer sequences to be removed, and no if you don’t. 
  * `THREADS (-d)` - This option lets you provide the number of threads available for use by the software. Default is 1. 
  * `NCBI_CLUSTERED_FILE (-b)` - This option helps obtain information about the database OTUs which clustered together.
  * `VERBOSE (-v)` - Delete intermediate files or not? 
  * `SPLIT_OTU (-o)` - Do you want to split reads into separate FASTA files based on the OTUs they match to? Answer in yes or no. 
  * `SPLIT_OTU_METHOD (-h)` - Do you want to split reads into separate OTUs before or after EE filtering? Works only if -o is yes. Answer in before or after. 


### Output files explained:  
![alt text](https://github.com/jpearl01/mcsmrt/blob/master/file_names.png "File names and description")

### Elaborate description for important files:
1) OTU table - File name: post_OTU_table.txt 
  * Column 1 - OTU ID - The are the names/id's of the centroid OTU's that were obtained during clustering. 
  * Column 2 to last column - Counts for number of reads in each sample. The number of columns depends on the number of samples which were processed. 

2) UTAX results file - File name: post_reads.utax
  * Column 1 - OTU ID with number of reads which mapped to that OTU during clustering.
  * Column 2 - Taxonomic assignments from domain to species level along with confidence for each level. 
  * Column 3 - Taxonomic assignment till maxinum confidence.
  * Column 4 - Strand at which the reference aligned with the read. 

3) Final results file - File name: post_final results - This is just a merge of OTU table, UTAX results and blast results. 

4) All reads info file - File name: pre_all_reads_info.txt
  * Column 1 - Read name
  * Column 2 - Basename, i.e., name of the sample to which that read belonged to
  * Column 3 - Number of CCS passes
  * Column 4 - Barcode used for this sample while pooling together for sequencing
  * Column 5 - Name of the sample to which that read belonged to, without barcode
  * Column 6 - Expected error before trimming primers
  * Column 7 - Expected error after trimming primers
  * Column 8 - Length of the read before trimming primers
  * Column 9 - Length of the read after trimming
  * Column 10 - Boolean for whether the read mapped to host genome or not
  * Column 11 - Boolean for whether the read had a forward primer
  * Column 12 - Boolean for whether the read had a reverse primer
  * Column 13 and 14 - Start and end coordinates for forward primer
  * Column 15 and 16 - Start and end coordinates for reverse primer
  * Column 17 - Orientation of the read originally. Note that after finding primers, reads are oriented to be uniform, i.e., 3'->5' direction. 
  * Column 18 - Note on the quality of primers that were matched
  * Column 19 - Number of primers which were present in that read
  * Column 20-33 - Taxonomic assignment for each level along with confidences. This is obtained by running utax on all the reads. 


### Tutorial:
This section is a walk through for how to use MCSMRT using sample data from BEI sequencing on PacBio. Before this, CCS and demultiplexing was run on the samples using the Reads of Insert (ROI) protocol from SMRT link. A filtering criteria of minimum number of 5 ccs passes and predicted minimum accuracy of 90 was used. Successful completion of the ROI protocol creates a results directory in a location/path which is based on how PacBio was configured. A directory similar to the one produced by ROI is in the tutorial folder in mcsmrt's GitHub page. You can then run each of these steps to learn how to use MCSMRT.  


1) cd to your home directory and clone the MCSMRT repository using  
   `$ cd ~`    
   `$ git clone git@github.com:jpearl01/mcsmrt.git`   
   This will create a directory called mcsmrt in your home directory.   
2) Create a new directory called "mcmsrt_tutorial" in your home directory by running  
   `$ mkdir ~/mcmsrt_tutorial`
   This is where all the analysis will be carried out and results files will be stored.  
3) Enter the mcmsrt_tutorial_bei folder using  
   `$ cd ~/mcmsrt_tutorial`
3) Download the tutorial data from:
   https://drive.google.com/open?id=1UJZBU3PhEVq8lUGcjPcs2s2LbqjsQctA
   and place all files in the mcmsrt_tutorial folder. Use the unzip command if you downloaded all the files at once. Make sure the files from the tutorial are not embedded in a folder, they should be present as individual files in the mcsmrt_tutorial directory. 
4) Expand the archive BEI_sample_datatar.gz. It will create a folder with the data files from the PacBio's ROI protocol. This data is from a single cell of a PacBio sequencing run, and contain 4 barcoded replicates of the BEI mock community. Uncompress using:  
   `$ tar -xzf BEI_sample_data.tar.gz`
   Successful completion should create a folder called data in your current working directory. This folder is structured exactly as any ROI protocol output file.  
5) In the sample_key.txt file, change the data_path column to the FULL PATH where your data folder is located, i.e., the folder that was created after uncompressing the tar.gz file. Do not use relative paths here. Again, this is the folder which looks like an ROI results folder.  
6) Run the get_fastqs.rb script to obtain FASTQ files with modified headers.  
   `$ ruby ~/mcsmrt/get_fastqs.rb -s sample_key.txt -o reads`
   On successful completion, a folder called reads should be created along with the FASTQ files.  
7) The last step is to run the mcsmrt.rb script which does most of the heavy lifting and produces an OTU table with taxonomies for each OTU. The command that should be run for this purpose is as follows  
   `$ ruby ~/mcsmrt/mcsmrt.rb -a -f reads/ -d num_of_threads_available -e 1 -s 5 -x 2000 -n 500 -c ~/mcsmrt/tutorial/rdp_gold.fa -t ~/mcsmrt/tutorial/16sMicrobial_ncbi_lineage_reference_database.udb -l ~/mcsmrt/tutorial/16sMicrobial_ncbi_lineage.fasta -g /path/to/human/genome/fasta -p ~/mcsmrt/tutorial/primers.fasta -b ~/mcsmrt/ncbi_clustered_table.tsv -v`
   With the `-d` option, provide the number of threads. With the `-g` option, provide the path to the complete human genome in FASTA format. The other input files required to run this script are provided in the tutorial folder.

### Example output files:
After successful completion of the command which runs the mcsmrt_v1.rb script, many output files are generated. The most imporatant/useful output files are pre_all_reads_info.txt and post_final_results.txt. As an example, results run through MCSMRT using data from BEI mock community is added in the folder called example_output_files. 


### Built with:  
ruby 2.2.1p85 (2015-02-26 revision 49769) [x86_64-linux]
