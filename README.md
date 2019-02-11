[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [](#lang-us) ![ruby in bioinformatics ftw](https://img.shields.io/badge/Language-ruby-steelblue.svg)


# MCSMRT
## Microbiome Classifier for SMRT PacBio data

### Citation:
>Earl, Joshua P., Nithin D. Adappa, Jaroslaw Krol, Archana S. Bhat, Sergey Balashov, Rachel L. Ehrlich, James N. Palmer, et al. 2018. “Species-Level Bacterial Community Profiling of the Healthy Sinonasal Microbiome Using Pacific Biosciences Sequencing of Full-Length 16S rRNA Genes.” Microbiome 6 (1): 190.

### Introduction:
MCSMRT is a tool to cluster PacBio FL16S amplicon microbiome sequences into Operational Taxonomic Units (OTU) and assign species level taxonomic classifications. Outputs include a table of read counts assigned to each OTU centroid sequence (per sample) with corresponding taxonomic lineage, and a a table of read specific metrics (e.g., CCS count, expected error, length, primer matching result, etc.). 

### Dependencies: 

1. [Ruby](https://www.ruby-lang.org/en/) v2.2.1 or greater 
   
   Fedora:
   ```
   $ sudo dnf install ruby
   ```

   Centos/RHEL
   ```
   $ sudo yum install ruby
   ```

2.  [Bundler Gem](https://bundler.io/): 
    ```
    $ gem install bundler
    ```

3.  Other Gems: 
    ```
    $ bundle
    ```

4. [BWA](https://sourceforge.net/projects/bio-bwa/files/) (install via site instructions).
  Add to path. e.g.:
    ```
    $ ln -s /path/to/bwa ~/bin/
    ```

5. [Sambamba](http://lomereiter.github.io/sambamba/)
   Precompiled binary. Add to path, e.g.:
   ```
   $ tar xvf sambamba* && ln -s /path/to/sambamba ~/bin/
   ```


6. [Usearch v8.1](http://www.drive5.com/usearch/download.html) free 32bit version (64 bit preferred, but not free). Make executable:
   ```
   $ chmod +x usearch8.1.*
   ```
   Add to path (executable *must* be called 'usearch'). e.g.:
   ```
   $ln -s /home/user/apps/usearch8.1.1861_i86linux32 ~/bin/usearch
   ```

8. [h5py](https://www.h5py.org/) (for ccs_passes.py):
   ```
   $ sudo dnf install h5py
   ```

### Installation:
`$ git clone git@github.com:jpearl01/mcsmrt.git`

### Data Prerequisites:
Data | File type | Description
--- | --- | ---
PCR Primer Sequences | [Fasta](https://en.wikipedia.org/wiki/FASTA_format) | Forward and reverse PCR primer sequences.
[Taxonomy Classification Database](https://drive.google.com/open?id=1UJZBU3PhEVq8lUGcjPcs2s2LbqjsQctA) | [UDB](https://www.drive5.com/usearch/manual/udb_files.html) | USEARCH formatted species level taxonomy classification database* Download the [64bit](https://www.dropbox.com/s/2squcpwg3lx19od/16S_NCBI_64bit_utax8.1.1861.udb?dl=0) or [32bit](https://www.dropbox.com/s/g5tjy7unezde1o3/16S_NCBI_32bit_utax8.1.1861.udb?dl=0) NCBI FL16S db version for usearch v8.1.861 (**Please Note:** you must use the same formatted db as the usearch executable you are using)
[Taxonomy Classification Database](https://drive.google.com/open?id=1UJZBU3PhEVq8lUGcjPcs2s2LbqjsQctA) | [Fasta](https://en.wikipedia.org/wiki/FASTA_format) | As above, except formatted in FASTA. [Download](https://www.dropbox.com/s/kje22s4gdad2fkp/16S_NCBI_utax_and_sintax_formatted.fasta.gz?dl=0). See [Lineanator](https://github.com/bhatarchanas/lineanator) for directions to create your own.  
Clustered Tax DB | [TSV](https://en.wikipedia.org/wiki/Tab-separated_values) | Table of cluster assignments that defines the number of closely related species to each entry in the database. [Download](ncbi_clustered_table.tsv) for NCBI FL16S db above
[RDP gold database](http://drive5.com/uchime/gold.fa) | [Fasta](https://en.wikipedia.org/wiki/FASTA_format) | Trusted sequences used to identify chimeras with uchime. [Download](http://drive5.com/uchime/gold.fa)
[Host Genome](https://useast.ensembl.org/Homo_sapiens/Info/Index) | [Fasta](https://en.wikipedia.org/wiki/FASTA_format) | Fasta of host genome file, indexed with BWA** [Human Genome Download](https://useast.ensembl.org/Homo_sapiens/Info/Index)

*[Linnaean](https://en.wikipedia.org/wiki/Linnaean_taxonomy) taxonomy to classify OTU sequences. Linked databases were created using a curated set of the [NCBI's](https://www.ncbi.nlm.nih.gov/) 16S BLAST and taxonomy databases. To create your own [see here](http://www.drive5.com/usearch/manual/cmd_makeudb_utax.html). A tool that can generate this database called Lineanator is available [here](https://github.com/bhatarchanas/lineanator). 

**This may or may not make sense, depending on the provenance of your samples. This file is for filtration of off-target host sequence reads.

## Tutorial and Links to Example Data 
[**Find Here**](tutorial.md)

## Getting Started With Your Own Data
### Method 1: SMRTportal/SMRTlink demultiplexed CCS fastq files (one file per sample):
---
  * Run get_fastqs.rb. This script requires fastq files from a Reads of Insert (ROI) protocol demultiplexed using SMRT portal/SMRT link, and the original output folders/data. This will insert CCS counts and the barcode name into the fastq headers. Modified fastq files are deposited in a user defined output folder. You will need a 'sample_info' file describing each samples name, and its corresponding output folder ([format described below](#sample_info-file-format)).

Command:    
```
  get_fastqs.rb 
  [-h] Help
  [-o OUTPUT_FOLDER_NAME] Folder for fastq file output. Fastq filenames become the sample_name from SAMPLE_INFO_FILE. 
  [-s SAMPLE_INFO_FILE] Tab separated file of all fastq files by sample.
```
 Example:

`$ ruby get_fastqs.rb -s samp_info.tsv -o mcsmrt_reads`

SAMPLE_INFO file format (TSV) Column detail:

Column | Description
--- | --- 
PB_jobid | Job ID assigned during demultiplexing (via SMRT portal).
data_path | Path to SMRT portal demultiplexed data.
forward_barcode | Forward primer name, (must match header in PRIMERS_DB file).
reverse_barcode | Reverse primer name, (must match header in PRIMERS_DB file).
sample_name | Sample/barcode label (added to fastq header as “barcodelabel”). Each unique forward/reverse barcode pair is given a different number that is added to the beginning of the “barcodelabel”. **Please Note**: Do not use characters that will break filenames, e.g. spaces, reserved system characters etc. 

Example formatted file:

PB_jobid | data_path | forward_barcode | reverse_barcode | sample_name
--- | --- | --- | --- | --- 
17094 | /data/pacbio/smrtanalysis_userdata/jobs/017/017094/ | 0001_Forw_bc | 0002_Rev_bc | A_22_gotaq
17094 | /data/pacbio/smrtanalysis_userdata/jobs/017/017094/ | 0003_Forw_bc | 0004_Rev_bc | A_hum_22_gotaq



### or

### Method 2: Demultiplex samples by barcode from raw pacbio output
---
  * Demultiplex using [this demultiplexing pipeline](https://github.com/rehrlich/ccs_smrt_pipe). This software produces FASTQ files formatted for the mcsmrt.rb script. This requires you have Smrtanalysis v2.3 installed

### Usage: 

 **Run mcsmrt.rb:**
    
   Either provide a directory with all fastq files for clustering, or provide a file with a list of file paths.

   Command parameters:
  
  ```
    mcsmrt.rb [-h]  
  [-a] / [-i LIST_OF_FILES_FOR_CLUSTERING] 
    [-f FOLDER_NAME] 
    [-m TRIMMING]  
  [-e EXPECTED_ERROR] 
    [-s CCS_COUNT] 
    [-x MAXIMUM_LENGTH] 
    [-n MINIMUM_LENGTH]  
  [-c UCHIME_DB] 
    [-t UTAX_DB] 
    [-l BLAST_DB] 
    [-g HOST_GENOME] 
    [-p PRIMERS]
  [-d THREADS] 
    [-b NCBI_CLUSTERED_FILE] 
    [-v VERBOSE]
  ``` 
   Example (directory method):


```
$ mcsmrt.rb -d 32 -a -f reads/ \
-c rdp_gold.fa \
-t 16S_NCBI_32bit_utax8.1.1861.udb \
-l 16S_NCBI_utax_and_sintax_formatted.fasta \
-g human_g1k_v37.fasta \
-p primers.fasta \
-b ncbi_clustered_table.tsv -v 
```
   Example (file method):

```
$ mcsmrt.rb -d 32 -i paths_to_fastq.tsv \
-c rdp_gold.fa \
-t 16S_NCBI_32bit_utax8.1.1861.udb \
-l 16S_NCBI_utax_and_sintax_formatted.fasta \
-g human_g1k_v37.fasta \
-p primers.fasta \
-b ncbi_clustered_table.tsv -v 
```


#### mcsmrt.rb Parameter Details:
Option | Description
--------- | ----------- 
**Mandatory** | 
`ALL (-a)` | Folder with all fastq files to use in clustering (assumes filename matches sample name)
`LIST_OF_FILES_FOR_CLUSTERING (-i)` | Specify sequence files with a file of fastq file paths, one file/location per row, tab separated
`FOLDER_NAME (-f)` | Project folder path. The input fastqs must be in a sub folder here.
`HOST_GENOME_DB (-g)` | Fasta formatted host reference genome file, indexed with BWA. Used to filter off-target host-mapping reads.
`PRIMERS_DB (-p)` | Fasta format primer file for primer matching and trimming.
`UCHIME_DB (-c)` | Database file with trusted high quality 16s sequences (fasta format). Default is rdp_gold database. Used to filter chimeric OTU sequences.
`UTAX_DB (-t)` | UDB format database file with lineage assigned to each 16S sequence. Custom database file was created with [Lineanator](https://github.com/bhatarchanas/lineanator) using the NCBI 16S BLAST and taxonomy databases. To learn more, [see this reference](http://www.drive5.com/usearch/manual/cmd_makeudb_utax.html).
`BLAST_DB (-l)` | The fasta file used to create the `UTAX_DB` file. Used to obtain strain name, alignment length, and percent identity for each OTU sequence.
`NCBI_CLUSTERED_FILE (-b)` | Table of database OTUs (dbOTU)
**Optional** |
`EXPECTED_ERROR (-e)` | (Default 1.0) Maximum [expected error](http://www.drive5.com/usearch/manual/expected_errors.html) (EE). Reads with higher EE are removed from further analysis.
`CCS_COUNT (-s)` | Default (5) – Minimum CCS count below which sequences are filtered out.
`MAXIMUM_LENGTH (-x)` | Default (2000) – Maximum length above which sequences are filtered out.
`MINIMUM_LENGTH (-n)` | Default (500) – Minimum length below which sequences are filtered out.
`TRIMMING (-m)` | Default (yes) – `yes` to trim primer sequences, and `no` otherwise.
`THREADS (-d)` | Number of threads to use. Default is 1.
`VERBOSE (-v)` | Keep all intermediate files
`SPLIT_OTU (-o)` | Default (no) - `yes` to split reads mapping to each OTU into separate multi-FASTA files. `no` - don't
`SPLIT_OTU_METHOD (-h)` | When `SPLIT_OTU` is `yes`, further define splitting OTU mapping reads before or after EE filtering. `before` or `after`. 


For a detailed description of the output files, please see [here](detailed_file_outputs.md)




### Built with:  
ruby 2.2.1p85 (2015-02-26 revision 49769) [x86_64-linux]
