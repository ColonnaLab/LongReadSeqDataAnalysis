[back to course home page ](../README.md)

## Variant Calling 
This page contain a short and adapted version of the Software Carpentry lesson [Variant Calling Workflow ](https://datacarpentry.github.io/wrangling-genomics/04-variant_calling.html). We will use it as notes to key concepts we will discuss during our lesson

This workflow demonstrates how to identify genetic variants from aligned sequencing data. Each step transforms the data from raw alignments to filtered variant calls ready for downstream analysis. The pipeline uses standard tools like BWA for alignment, SAMtools for processing, and bcftools for variant calling.


#### getting the reference genome 
Before calling variants, we need to align our reads to a reference genome. we download the reference genome for E. coli REL606

Within `seq-analysis` make a dir called refgenome and navigate there 

```diff
user1@vm-corso-colonna:~/seq-analysis$ mkdir refgenome
user1@vm-corso-colonna:~/seq-analysis$ cd refgenome/
user1@vm-corso-colonna:~/seq-analysis/refgenome$ 
```

Download, gunzip, index 

```diff 

+ COMMAND  curl - transfer a URL

user1@vm-corso-colonna:~/seq-analysis/refgenome$ curl -L -o ecoli_rel606.fasta.gz \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz

+ COMMAND gunzip 
 user1@vm-corso-colonna:~/seq-analysis/refgenome$ gunzip ecoli_rel606.fasta.gz

```

#### alignment with `bwa`

**[BWA](http://bio-bwa.sourceforge.net/)** (Burrows-Wheeler Aligner) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome, using the Burrows-Wheeler transform algorithm.

The reference genome must be indexed before alignment:

```bash
user1@vm-corso-colonna:~/seq-analysis/variant-calling$ bwa index ../data/ref_genome/ecoli_rel606.fasta
user1@vm-corso-colonna:~/seq-analysis/refgenome$ bwa index ecoli_rel606.fasta 
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.18 seconds elapse.
[bwa_index] Update BWT... 0.04 sec
[bwa_index] Pack forward-only FASTA... 0.03 sec
[bwa_index] Construct SA from BWT and Occ... 0.46 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index ecoli_rel606.fasta
[main] Real time: 2.060 sec; CPU: 1.755 sec
user1@vm-corso-colonna:~/seq-analysis/refgenome$ ls 
ecoli_rel606.fasta      ecoli_rel606.fasta.ann  ecoli_rel606.fasta.pac
ecoli_rel606.fasta.amb  ecoli_rel606.fasta.bwt  ecoli_rel606.fasta.sa

```

```diff
! Why do we index? Indexing creates a searchable data structure that speeds up alignment
```

Go to `~/seq-analysis` Align paired-end reads using BWA-MEM:

```bash
user1@vm-corso-colonna:~/seq-analysis$ bwa mem \
    refgenome/ecoli_rel606.fasta \
    trimmed_all/SRR2584863_1.trimmed.fastq 
    trimmed_all/SRR2584863_2.trimmed.fastq \
    > results/sam/SRR2584863.aligned.sam 
```

```diff
+ BWA MEM is optimized for reads 70bp-1Mbp
+ Outputs SAM format to stdout (we redirect to file with >)
+ Uses paired-end mode when given two input files
```

##### preparing the filesystem for results 

mkdir res 

```
user1@vm-corso-colonna:~/seq-analysis$ mkdir -p results/sam results/bam results/bcf results/vcf
user1@vm-corso-colonna:~/seq-analysis$ ls 
fastqc-res  refgenome  results  trimmed  trimmed_all
user1@vm-corso-colonna:~/seq-analysis$ ls results/
bam  bcf  sam  vcf
user1@vm-corso-colonna:~/seq-analysis$ 
```



