require 'bio'
require 'trollop'

##### Input 
opts = Trollop::options do
	opt :eevalue, "Expected error at which you want filtering to take place", :type => :string, :short => "-e"
	opt :otufastafile, "Path for file with all the OTU sequences in FASTA format", :type => :string, :short => "-o"
	opt :spoofeddbfile, "Path for binary DB file for the spoofed data", :type => :string, :short => "-d"
	opt :lineagefastafile, "Path of FASTA file with lineage info for the ublast command", :type => :string, :short => "-l"
end 

##### Assigning variables to the input
ee = opts[:eevalue]
otu_fasta_file = opts[:otufastafile]
spoofed_db_file = opts[:spoofeddbfile]
lineage_fasta_file = opts[:lineagefastafile]

# For assigning taxonomy, get the udb file using makeudb_utax command, use the spoofed db...
`usearch -utax #{otu_fasta_file} -db #{spoofed_db_file} -utaxout all_ee#{ee}_spoofed_reads.utax -utax_cutoff 0.8 -strand both`

# Running the command to give a report of counts
`ruby get_report_spoofed.rb -s all_ee#{ee}`

# Running blast on the OTUs                                                                                                                            
`usearch -ublast #{otu_fasta_file} -db #{lineage_fasta_file} -top_hit_only -id 0.9 -blast6out all_ee#{ee}_spoofed_blast.txt -strand both -evalue 0.01 -threads 15 -accel 0.3`

# Running the script which maps between the blast file and the merged OTUs file                                                                         
`ruby map_blast_with_otus.rb -u all_ee#{ee}_spoofed_table_utax_map.txt -b all_ee#{ee}_spoofed_blast.txt -o all_ee#{ee}_spoofed_otu_blast_formatted.txt`

# Running the script which dismantles the taxonomic assignments field
`ruby parse_final_header.rb -h all_ee#{ee}_spoofed_otu_blast_formatted.txt -o all_ee#{ee}_spoofed_otu_blast_formatted2.txt`  
