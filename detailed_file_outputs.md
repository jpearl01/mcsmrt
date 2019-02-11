### Output files:  
![alt text](images/file_names.png "File names and description")

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
  * Column 20-33 - Taxonomic assignment for each level along with confidences.  
