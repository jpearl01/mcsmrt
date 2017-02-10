##MCSMRT

Author - Archana S Bhat

###Introduction:

MCSMRT is a tool that can be used for microbiome analysis of PacBio data. Some of the information that can be obtained from this tool includes clustering into OTUs and counts for the number of reads belonging to each OTU, taxanomic lineage associated with each OTU. Another file that is produced as an output consists of general information regarding each and every read, for example, CCS count, expected error, length, primer matching result, etc. 


###Installation and Dependencies: 

1) Ruby gems such as bio and troloop are required for Lineanator to function. If they are not already installed, use the "gem install {name_of_the_gem}" command to install these gems.

2) Download and install usearch v8.1.
   Create a soft link pointing towards this version of usearch and name it "usearch". The soft link name **HAS TO BE** usearch, this is **VERY IMPORTANT**.


###Usage: 

There are 2 methods/routes you can choose from, in order to run your PacBio data through the microbiome classifier. 

Method 1: In case you have already demultiplexed your data
* Step 1: Run get_fastqs.rb in order to obtain the CCS and barcode information along with the demultiplexed FASTQ sequences in separate files within a folder name given by you. <br /><br /> ruby get_fastqs.rb [-h] [-s SAMPLE_INFO_FILE] [-o OUTPUT_FOLDER_NAME] 

* Step 2: Run mcsmrt.rb in order to obtain detailed information about each read and get the clustering information about the OTUs created. <br />
This step in-turn can be used in 2 ways, one way is to speficy that all the files in a given directory can be used for clustering and another way is to provide a file with a list of file names which you want to be clustered together. <br /><br /> ruby mcsmrt.rb [-h] [-i all_bc_reads.fastq] [-e 1] [-s 5] [-x 2000] [-n 500] [-c ../rdp_gold.fa] [-t ../lineanator/16sMicrobial_ncbi_lineage_reference_database.udb] [-l ../lineanator/16sMicrobial_ncbi_lineage.fasta] [-g ../human_g1k_v37.fasta] [-p ../primers.fasta]
