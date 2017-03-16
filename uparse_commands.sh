# For dereplication, -sizeout is the number of reads that are collapsed into one read, minseqlength option can be used to specify the minimum sequence length to be included in the output. 
usearch -derep_fulllength $1.fq -fastaout $1_dereplicated_s1.fa -sizeout 

# For clustering the OTU's, -otus option stands for the name of the output file in fasta format, -uparseout stands for the name of the file from uparse, -sizein and -sizeout says that size annotation is present in the input file and is required in the output file.
# id - http://drive5.com/usearch/manual/opt_id.html, specifies the minimum sequence identity of a hit.
# -otu_radius_pct 1.35 (cuz a paper josh read said that species can be differentiated at 98.65% cutoffs)
usearch -cluster_otus $1_dereplicated_s1.fa -otus $1_OTU_s2.fa -uparseout $1_uparse_s2.up -relabel OTU_ -sizein -sizeout 

# For filtering the chimeric reads, using the gold database, -nonchimeras needs the name of the file without the chimeric reads                              
usearch -uchime_ref $1_OTU_s2.fa -db $2 -strand plus -nonchimeras $1_OTU_nonchimeras_s3.fa -chimeras $1_OTU_chimeras_s3.fa -uchimealns $1_OTU_alignment_s3.aln -uchimeout $1_OTU_uchime_s3.txt

# For mapping the reads to the OTUS's, -id 0.97 says that we want 97% similarity between the reads and the reference sequence in the OTU list, in order for mapping to occur
usearch -usearch_global $4 -db $1_OTU_nonchimeras_s3.fa -strand plus -id 0.97 -uc $1_readmap_s4.uc -otutabout $1_OTU_table_s4.txt -notmatched $1_unmapped_userach_global.fa

# For assigning taxonomy, get the udb file using makeudb_utax command 
usearch -utax $1_OTU_nonchimeras_s3.fa -db $3 -utaxout $1_reads_s5.utax -utax_cutoff 0.8 -strand both

