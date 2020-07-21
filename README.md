[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [](#lang-us) ![ruby in bioinformatics ftw](https://img.shields.io/badge/Language-ruby-steelblue.svg)


# MCSMRT
## Microbiome Classifier for SMRT PacBio data

### Citation:
>Earl, Joshua P., Nithin D. Adappa, Jaroslaw Krol, Archana S. Bhat, Sergey Balashov, Rachel L. Ehrlich, James N. Palmer, et al. 2018. “Species-Level Bacterial Community Profiling of the Healthy Sinonasal Microbiome Using Pacific Biosciences Sequencing of Full-Length 16S rRNA Genes.” Microbiome 6 (1): 190.

### Introduction:
MCSMRT is a tool to cluster PacBio FL16S amplicon microbiome sequences into Operational Taxonomic Units (OTU) and assign species level taxonomic classifications. Outputs include a table of read counts assigned to each OTU centroid sequence (per sample) with corresponding taxonomic lineage, and a a table of read specific metrics (e.g., CCS count, expected error, length, primer matching result, etc.). 

### Installation:
```bash
  $ git clone git@github.com:jpearl01/mcsmrt.git
  $ cd mcsmrt
```

### Dependencies: 
(execute in the mcsmrt/ directory)

1. [Ruby](https://www.ruby-lang.org/en/) v2.2.1 or greater 
   
  The simplest way is to download and install the [Ruby Version Manager (rvm)](https://rvm.io/rvm/install) 

  then:

   ```bash
   $ rvm install 2.2.1
   $ rvm use 2.2.1
   ```

2.  [Bundler Gem](https://bundler.io/): 
    ```bash
    $ gem install bundler
    ```

3.  All Required Ruby Gems: 
    ```bash
    $ bundle
    ```

4. [BWA](https://sourceforge.net/projects/bio-bwa/files/) (install via site instructions).
  Add to path. e.g.:
    ```bash
    $ ln -s /path/to/bwa ~/bin/
    ```

5. [Sambamba](http://lomereiter.github.io/sambamba/)
   Precompiled binary. Add to path, e.g.:
   ```bash
   $ tar xvf sambamba* && chmod +x /path/to/sambamba && ln -s /path/to/sambamba ~/bin/sambamba
   ```


6. [Usearch v8.1](http://www.drive5.com/usearch/download.html) free 32bit version (64 bit preferred, but not free). Make executable:
   ```bash
   $ chmod +x usearch8.1.*
   ```
   Add to path (executable *must* be called 'usearch'). e.g.:
   ```bash
   $ ln -s /home/user/apps/usearch8.1.1861_i86linux32 ~/bin/usearch
   ```

8. [h5py](https://www.h5py.org/) (for ccs_passes.py):
   ```bash
   $ sudo dnf install h5py
   ```

 9. Python 2.7 (only required if using older RSII data for ccs_passes.py)
    [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) is a good way to manage python versions



### Data Prerequisites:
**UPDATE:** You may download all required files [here](https://www.dropbox.com/sh/xodhq7pek6rfkaa/AABKUuUJqeEE5pmWg9OwPo6oa?dl=0)
Data | File type | Description
--- | --- | ---
PCR Primer Sequences | [Fasta](https://en.wikipedia.org/wiki/FASTA_format) | Forward and reverse PCR primer sequences (currently must be a fasta file with two records named 'forward' and 'reverse').
[Taxonomy Classification Database](https://drive.google.com/open?id=1UJZBU3PhEVq8lUGcjPcs2s2LbqjsQctA) | [UDB](https://www.drive5.com/usearch/manual/udb_files.html) | USEARCH formatted species level taxonomy classification database* Download the [64bit](https://www.dropbox.com/s/2squcpwg3lx19od/16S_NCBI_64bit_utax8.1.1861.udb?dl=0) or [32bit](https://www.dropbox.com/s/g5tjy7unezde1o3/16S_NCBI_32bit_utax8.1.1861.udb?dl=0) NCBI FL16S db version for usearch v8.1.861 (**Please Note:** you must use the same formatted db as the usearch executable you are using the free version is **32bit**)
Clustered Tax DB | [TSV](https://en.wikipedia.org/wiki/Tab-separated_values) | Table of cluster assignments that defines the number of closely related species to each entry in the database. [Download](data/ncbi_clustered_table.tsv) for NCBI FL16S db above
[RDP gold database](http://drive5.com/uchime/gold.fa) | [Fasta](https://en.wikipedia.org/wiki/FASTA_format) | Trusted sequences used to identify chimeras with uchime. This file is called 'gold.fa' at this download link, and is equivalent to 'rdp_gold.fa' used in the example commands below [Download](http://drive5.com/uchime/gold.fa)
[Host Genome](https://useast.ensembl.org/Homo_sapiens/Info/Index) | [Fasta](https://en.wikipedia.org/wiki/FASTA_format) | Fasta of host genome file, MUST index with BWA after download** [Human Genome Download](https://useast.ensembl.org/Homo_sapiens/Info/Index)

*[Linnaean](https://en.wikipedia.org/wiki/Linnaean_taxonomy) taxonomy to classify OTU sequences. Linked databases were created using a curated set of the [NCBI's](https://www.ncbi.nlm.nih.gov/) 16S BLAST and taxonomy databases. To create your own [see here](http://www.drive5.com/usearch/manual/cmd_makeudb_utax.html). A tool that can generate this database called Lineanator is available [here](https://github.com/bhatarchanas/lineanator). 

**Using a host genome may or may not make sense, depending on the provenance of your samples. This file is for filtration of off-target host sequence reads. If your samples come from an organism other than Human, you would want do download the fasta reference genome of your other host.  If your samples come from the environment, or other non-host associated origin, you would not need this. Please see [indexing with BWA here](http://bio-bwa.sourceforge.net/bwa.shtml) to prepare your fasta reference for use with MCSMRT.  You probably want to use a relatively recent version of your host genome, but otherwise the version of the genome doesn't matter.

## Tutorial and Links to Example Data 
[**Find Here**](tutorial.md)

## Getting Started With Your Own Data
### Method 1: SMRTportal/SMRTlink demultiplexed CCS fastq files (one file per sample):
---
  * Run get_fastqs.rb. This script requires fastq files from a Reads of Insert (ROI) protocol demultiplexed using SMRT portal/SMRT link, and the original output folders/data. This will insert CCS counts and the barcode name into the fastq headers. Modified fastq files are deposited in a user defined output folder. You will need a 'sample_info' file describing each samples name, and its corresponding output folder (format described below).

Command:    

get_fastqs.rb Parameters:

Option | Description
--------- | ----------- 
`-h` | Help
`-o OUTPUT_FOLDER_NAME` | Folder for fastq file output. Fastq filenames become the sample_name from SAMPLE_INFO_FILE. 
`-s SAMPLE_INFO_FILE` | Tab separated file of all fastq files by sample.

 Example:

`$ ruby get_fastqs.rb -s samp_info.tsv -o mcsmrt_reads`

SAMPLE_INFO file format (TSV) Column detail:

Column_ID | Description
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


   Example (directory method):


```
$ mcsmrt.rb -d 32 -f reads/ \
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
`--fastqFolder\ (-f)` | Folder path of all fastq files.
`--samplelist\ (-i)` | File with a subset of tab separated sample/file_paths (one per row)
`--hostDB\ (-g)` | Fasta formatted host reference genome file, indexed with BWA. Used to filter off-target host-mapping reads.
`--primerfile\ (-p)` | Fasta format of all primer sequences. Needed for primer matching/trimming.
`--uchimedbfile\ (-c)` | Fasta formatted database file with trusted 16s sequences.  Used to filter chimeric OTU sequences.
`--utaxdbfile\ (-t)` | UDB format database file with lineage assigned to each 16S sequence. Custom database file was created with [Lineanator](https://github.com/bhatarchanas/lineanator) using the NCBI 16S BLAST and taxonomy databases. To learn more, [see this reference](http://www.drive5.com/usearch/manual/cmd_makeudb_utax.html).
`--clustereddb\ (-b)` | Table of database OTUs (dbOTU). Used to identify classifications with many closely related entries in the taxonomy database.
**Optional** |
`--eevalue\ (-e)` | (Default 1.0) Maximum [expected error](http://www.drive5.com/usearch/manual/expected_errors.html) (EE). Reads greater EE are removed.
`--ccsvalue\ (-s)` | Default (5) – Minimum CCS count below which sequences are filtered out.
`--lineagefastafile\ (-l)` | The fasta file used to create the `utaxdb` file. Used to obtain strain name, alignment length, and percent identity for each OTU sequence. **NO LONGER REQUIRED, DEPRECATED** 
`--lengthmax\ (-x)` | Default (2000) – Maximum read length.
`--lengthmin\ (-n)` | Default (500) – Minimum read length.
`--trimming\ (-m)` | Default (yes) – `yes` to trim primer sequences, and `no` otherwise.
`--threads\ (-d)` | Default (1) - Number of threads to use.
`--verbose\ (-v)` | Keep all intermediate files
`--splitotu\ (-o)` | Default (no) - `yes` to split reads mapping to each OTU into separate multi-FASTA files. `no` - don't
`--splitotumethod\ (-j)` | When `splitotu` is `yes`, further define splitting OTU mapping reads before or after EE filtering. `before` or `after`. 


For a detailed description of the output files, please see [here](detailed_file_outputs.md)




### Built with:  
ruby 2.2.1p85 (2015-02-26 revision 49769) [x86_64-linux]
