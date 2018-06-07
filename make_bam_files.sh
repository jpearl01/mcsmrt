bwa index first_16s_copy.fa
bwa mem -t 16 first_16s_copy.fa ccs2_First\&Sixth_subset.fa > map_to_first_16s.sam
sambamba view -S -f bam map_to_first\&sixth_16s.sam -o map_to_first\&sixth_16s.bam
samtools mpileup -uf first_16s_copy.fa map_to_first\&sixth_16s.bam | bcftools call -mv -o map_to_first\&sixth_16s.vcf