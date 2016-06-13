# For filtering the low quality reads, -fastq_qmax 50 says that 50 is the maximum quality possible, -fastq_maxee_rate 1.0 says that the reads which have expected error more than 1.0 are to be removed. 
usearch -fastq_filter $1_oriented_primer_trimmed_final.fastq -fastqout $1_oriented_filtered.fastq -fastq_maxee $2 -fastq_eeout -sample all

# For dereplication, -sizeout is the number of reads that are collapsed into one read, minseqlength option can be used to specify the minimum sequence length to be included in the output. 
usearch -derep_fulllength $1_oriented_filtered.fastq -fastaout $1_oriented_filtered_derep.fasta -sizeout 

# For clustering the OTU's, -otus option stands for the name of the output file in fasta format, -uparseout stands for the name of the file from uparse, -sizein and -sizeout says that size annotation is present in the input file and is required in the output file.
# id - http://drive5.com/usearch/manual/opt_id.html, specifies the minimum sequence identity of a hit.
# -otu_radius_pct 1.35 (cuz a paper josh read said that species can be differentiated at 98.65% cutoffs)
usearch -cluster_otus $1_oriented_filtered_derep.fasta -otus $1_OTU.fasta -uparseout $1_uparse_out.up -relabel OTU_ -sizein -sizeout 

# For filtering the chimeric reads, using the gold database, -nonchimeras needs the name of the file without the chimeric reads                              
usearch -uchime_ref $1_OTU.fasta -db $3 -strand plus -nonchimeras $1_OTU_nonchimera.fasta -chimeras $1_OTU_chimeras.fasta -uchimealns $1_OTU_alignment.aln -uchimeout $1_OTU_uchime.txt

# For mapping the reads to the OTUS's, -id 0.97 says that we want 97% similarity between the reads and the reference sequence in the OTU list, in order for mapping to occur
usearch -usearch_global $1_oriented_primer_trimmed_final.fastq -db $1_OTU.fasta -strand plus -id 0.97 -uc $1_readmap.uc -otutabout $1_Table.txt

# For assigning taxonomy, get the udb file using makeudb_utax command 
usearch -utax $1_OTU.fasta -db $4 -utaxout $1_reads.utax -utax_cutoff 0.8 -strand both

