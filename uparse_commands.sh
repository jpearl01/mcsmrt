#!/usr/bin/env bash
# For dereplication, -sizeout is the number of reads that are collapsed into one read, minseqlength option can be used to specify the minimum sequence length to be included in the output. 
usearch -derep_fulllength $1.fq -fastaout post_dereplicated.fa -sizeout 

# For clustering the OTU's, -otus option stands for the name of the output file in fasta format, -uparseout stands for the name of the file from uparse, -sizein and -sizeout says that size annotation is present in the input file and is required in the output file.
# id - http://drive5.com/usearch/manual/opt_id.html, specifies the minimum sequence identity of a hit.
# -otu_radius_pct 1.35 (cuz a paper josh read said that species can be differentiated at 98.65% cutoffs)
usearch -cluster_otus post_dereplicated.fa -otus post_OTU.fa -uparseout post_uparse.up -relabel OTU_ -sizein -sizeout 

# For filtering the chimeric reads, using the gold database, -nonchimeras needs the name of the file without the chimeric reads                              
usearch -uchime_ref post_OTU.fa -db $2 -strand plus -nonchimeras post_OTU_nonchimeras.fa -chimeras post_OTU_chimeras.fa -uchimealns post_OTU_alignment.aln -uchimeout post_OTU_uchime_output.tsv

# For mapping the reads to the OTUS's, -id 0.97 says that we want 97% similarity between the reads and the reference sequence in the OTU list, in order for mapping to occur
usearch -usearch_global $4 -db post_OTU_nonchimeras.fa -strand plus -id 0.97 -uc post_readmap.uc -otutabout post_OTU_table.tsv -notmatched post_unmapped_userach_global.fa -userout post_usearch_glob_results.tsv -userfields query+target

# For assigning taxonomy, get the udb file using makeudb_utax command 
usearch -utax post_OTU_nonchimeras.fa -db $3 -utaxout post_reads.utax -utax_cutoff 0.8 -strand both


