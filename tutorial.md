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
3) Download the tutorial data from [here](https://drive.google.com/open?id=1UJZBU3PhEVq8lUGcjPcs2s2LbqjsQctA)
   and place all files in the mcmsrt_tutorial folder.
4) Expand the archive BEI_sample_datatar.gz which will create a folder with the data files from the PacBio's ROI protocol. These data are from a single cell of a PacBio sequencing run, and contain 4 barcoded replicates of the BEI mock community. Uncompress via:  
   ```
   $ tar -xzf BEI_sample_data.tar.gz
   ```
   Successful completion should create a folder called data in your current working directory. This folder is structured exactly as any [SMRT portal/SMRT link](https://www.pacb.com/products-and-services/analytical-software/smrt-analysis/) ROI protocol output.  
5) In the sample_key.txt file, you must change the data_path column to the FULL PATH to the newly created data folder. Do not use relative paths here. Again, this is the folder which looks like an ROI results folder.  
   and place all files in the mcmsrt_tutorial folder. Use the unzip command if you downloaded all the files at once. Make sure the files from the tutorial are not embedded in a folder, they should be present as individual files in the mcsmrt_tutorial directory. 
4) Expand the archive BEI_sample_data.tar.gz. It will create a folder with the data files from the PacBio's ROI protocol. This data is from a single cell of a PacBio sequencing run, and contain 4 barcoded replicates of the BEI mock community. Uncompress using:  
   `$ tar -xzf BEI_sample_data.tar.gz`
   Successful completion should create a folder called data in your current working directory. This folder is structured exactly as any ROI protocol output file.  
5) In the sample_key.txt file, change the data_path column to the FULL PATH where your data folder is located, i.e., the folder that was created after uncompressing the tar.gz file. Do not use relative paths here. Again, this is the folder which looks like an ROI results folder.  

6) Run the get_fastqs.rb script to obtain FASTQ files with modified headers.  
   
   ```
   $ ruby ~/mcsmrt/get_fastqs.rb -s sample_key.txt -o reads
   ```

   On successful completion, a folder called reads should be created along with the FASTQ files.  
7) The last step is to run the mcsmrt.rb script which does most of the heavy lifting and produces an OTU table with taxonomies for each OTU. The command that should be run for this purpose is as follows  

   ```
   $ ruby ~/mcsmrt/mcsmrt_v1.rb -a -f reads/ -d num_of_threads_available -e 1 -s 5 -x 2000 -n 500 -c ~/mcsmrt/tutorial/rdp_gold.fa -t ~/mcsmrt/tutorial/16sMicrobial_ncbi_lineage_reference_database.udb -l ~/mcsmrt/tutorial/16sMicrobial_ncbi_lineage.fasta -g /path/to/human/genome/fasta -p ~/mcsmrt/tutorial/primers.fasta -b ~/mcsmrt/ncbi_clustered_table.tsv -v
   ```
   
   `$ ruby ~/mcsmrt/mcsmrt.rb -a -f reads/ -d num_of_threads_available -e 1 -s 5 -x 2000 -n 500 -c ~/mcsmrt/tutorial/rdp_gold.fa -t ~/mcsmrt/tutorial/16sMicrobial_ncbi_lineage_reference_database.udb -l ~/mcsmrt/tutorial/16sMicrobial_ncbi_lineage.fasta -g /path/to/human/genome/fasta -p ~/mcsmrt/tutorial/primers.fasta -b ~/mcsmrt/ncbi_clustered_table.tsv -v`

   With the `-d` option, provide the number of threads. With the `-g` option, provide the path to the complete human genome in FASTA format. The other input files required to run this script are provided in the tutorial folder.

### Example output files:
After successful completion of the command which runs the mcsmrt_v1.rb script, many output files are generated. The most imporatant/useful output files are pre_all_reads_info.txt and post_final_results.txt. As an example, results run through MCSMRT using data from BEI mock community is added in the folder called example_output_files. 

