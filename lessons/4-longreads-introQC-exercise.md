### Data set 
## Dataset description

In this practical, we will analyse public Oxford Nanopore long-read sequencing data from the methicillin-resistant bacterium *Staphylococcus aureus* (MRSA), strain KUN1163. MRSA is a drug-resistant form of *S. aureus*. The dataset is suitable for a bacterial genomics workflow because the genome is relatively small and has a high sequencing depth.

The nanopore dataset is available through the International Nucleotide Sequence Database Collaboration under run accession **DRR187567**, within study **DRA008776**. The reads were generated on an Oxford Nanopore MinION instrument using an R9 flow cell. The raw dataset contains **91,288 reads** comprising approximately **622 Mb** of sequence data. Since the KUN1163 chromosome is approximately **2.91 Mb**, this corresponds to about **214× long-read coverage**.

For reference-based analysis, we will use the published complete genome of the same isolate:

* Chromosome: **AP020324.1**; 2,914,567 bp
* Plasmid p01KUN1163: **AP020325.1**; 30,220 bp

The original publication produced the completed genome using both Oxford Nanopore long reads and Illumina short reads. In this class, we will use the nanopore data to perform read-quality control, inspect read-length and quality distributions, filter reads, align long reads to the reference chromosome and plasmid, assess alignment and coverage statistics, and perform variant calling.

**Dataset citation**

Hikichi M, Nagao M, Murase K, Aikawa C, Nozawa T, Yoshida A, Kikuchi T, Nakagawa I. 2019. Complete Genome Sequences of Eight Methicillin-Resistant *Staphylococcus aureus* Strains Isolated from Patients in Japan. *Microbiology Resource Announcements* 8. https://doi.org/10.1128/mra.01212-19

###  Downloading the sequencing data from the SRA using the SRA Toolkit

Install the NCBI SRA Toolkit, then create a project directory and download run **DRR187567**.

```bash
mkdir /data-longreads
cd /data-longreads


# Download the SRA run
prefetch DRR187567

# Convert the SRA archive to FASTQ format
fasterq-dump DRR187567 \
  --threads 4 \
  --outdir fastq \
  --temp tmp

# Compress the FASTQ file to save disk space
gzip fastq/DRR187567.fastq
```

The resulting read file will be:

```text
fastq/DRR187567.fastq.gz
```

This is a single-end Oxford Nanopore FASTQ file.


### Download the reference genome

Download both the chromosome and plasmid reference sequences:

```bash
curl -L --fail \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP020324.1,AP020325.1&rettype=fasta&retmode=text" \
  > reference/KUN1163_reference.fasta

samtools faidx reference/KUN1163_reference.fasta
```

The file `reference/KUN1163_reference.fasta` should be used for long-read alignment and coverage analysis.

### Quality Control with Nanoplot 

this is where the data lives 

```bash 
user1@vm-corso-colonna:~$ ls data-longreads
DRR187567  fastq  reference
user1@vm-corso-colonna:~$ ls data-longreads/fastq/
DRR187567.fastq.gz
user1@vm-corso-colonna:~$ ls data-longreads/reference/
KUN1163_reference.fasta  KUN1163_reference.fasta.fai
user1@vm-corso-colonna:~$ 
```

#### Practical: summarize read length, quality, and yield

- Goal of the practical:
  - Start from the Oxford Nanopore FASTQ file
  - Count the number of reads
  - Estimate total yield in bases
  - Summarize read length statistics
  - Summarize read quality statistics
  - Generate an HTML QC report with plots

```diff
! EXERCISE: Move to your home folder
! Check that the long reads data and reference files are present
```

```bash
user1@vm-corso-colonna:~$ cd
user1@vm-corso-colonna:~$ ls data-longreads
user1@vm-corso-colonna:~$ ls data-longreads/fastq/
user1@vm-corso-colonna:~$ ls data-longreads/reference/
```

- First inspect the input FASTQ file:
  - Check the file size
  - Confirm that the file is compressed
  - Confirm that it is a FASTQ file

```diff
+ COMMAND: ls -lh - show file sizes in human-readable format
+ COMMAND: file - guess file type
```

```bash
user1@vm-corso-colonna:~$ ls -lh data-longreads/fastq/DRR187567.fastq.gz
user1@vm-corso-colonna:~$ file data-longreads/fastq/DRR187567.fastq.gz
```

- Use `seqkit stats` for a first command-line summary:
  - Number of reads
  - Total number of bases
  - Shortest and longest read
  - Average read length
  - N50 read length

```diff
+ COMMAND: seqkit stats - summarize FASTA/FASTQ files
```

```bash
user1@vm-corso-colonna:~$ seqkit stats data-longreads/fastq/DRR187567.fastq.gz
```

```diff
! EXERCISE: From the seqkit output, write down:
! 1. number of reads
! 2. total yield in bases
! 3. minimum read length
! 4. maximum read length
! 5. average read length
! 6. N50
```

- Create a dedicated output folder for the NanoPlot report:

```bash
user1@vm-corso-colonna:~$ mkdir -p qc_nanoplot
```

- Use `NanoPlot` to generate plots and an HTML report:
  - Read length distribution
  - Quality score distribution
  - Yield summary
  - Length versus quality plot

```diff
+ COMMAND: NanoPlot - make QC plots for long reads data
+ --fastq: input FASTQ file
+ --outdir: output directory
+ --prefix: prefix for output files
```

```bash
user1@vm-corso-colonna:~$ NanoPlot \
  --fastq data-longreads/fastq/DRR187567.fastq.gz \
  --outdir data-longreads/qc_nanoplot \
  --prefix DRR187567_
```

- Inspect the output folder:
  - HTML report
  - Summary statistics
  - Plots in image format
  - Log files

```bash
user1@vm-corso-colonna:~$ ls -lh data-longreads/qc_nanoplot/
```

```diff
! EXERCISE: Open the NanoPlot HTML report
! Identify:
! 1. the read length distribution
! 2. the quality distribution
! 3. the total yield
! 4. whether there are many very short reads
! 5. whether long reads also have acceptable quality
```

- Questions to answer after the practical:
  - Is the dataset dominated by short, medium, or very long reads?
  - Are the longest reads also high quality?
  - Would you remove very short reads before alignment?
  - Would you remove low-quality reads before alignment?
  - Is the total yield sufficient for the expected genome size?

```diff
+ Remember: filtering decisions depend on the downstream analysis
+ For assembly, very long reads may be useful even if their quality is lower
+ For variant calling, quality and coverage are both important
```
