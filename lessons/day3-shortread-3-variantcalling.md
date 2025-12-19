[back to course home page ](../README.md)

## Variant Calling 
This page contain a short and adapted version of the Software Carpentry lesson [Variant Calling Workflow ](https://datacarpentry.github.io/wrangling-genomics/04-variant_calling.html). We will use it as notes to key concepts we will discuss during our lesson

This workflow demonstrates how to identify genetic variants from aligned sequencing data. Each step transforms the data from raw alignments to filtered variant calls ready for downstream analysis. The pipeline uses standard tools like BWA for alignment, SAMtools for processing, and freebayes for variant calling. Here is an overview of the vomplete workflow: 

```diff
+ Reads → BWA → SAM → SAMtools → BAM → FreeBayes → VCF → BCFtools → Filtered VCF
```

#### What Does It Mean to "Align" DNA Sequences?
The first step in variant calling is  the sequence alignment. **Aligning** means finding where a short piece of DNA fits best within a much longer DNA sequence - like finding where a sentence fragment belongs in a book. The computer tries millions of positions to find the best match, accounting for small differences that naturally occur between individuals.

```diff
Book (Reference):    "The quick brown fox jumps over the lazy dog"
                                ↑_____________↑
Sentence fragment:             "brown fox jumps"
(Your DNA read)

In DNA Terms

Reference DNA:  ...ATCGATCGATCGATCG...
                   ||||||||||||||||
Your DNA read:     ATCGATCGATCGATCG

Differences between reference and query: 
| **Mismatches** | Ref: `ATCG`  Read: `ATGG` |
| **Insertions** | Ref: `AT_CG` Read: `ATACG`|
| **Deletions**  | Ref: `ATCG`  Read: `AT_G` |

+ Key Concepts
+ 1. Reference Genome: it is the complete genomic sequence (template) and it is indexed for fast searching
+ 2. Query Read: it is a short (typically 50-300 bp) or long (typically 10-30 kbp) DNA fragment from sequencing Typically 50-300 bp

+ Short Reads:    100-300 bp       |-----|
+ Long Reads:     10,000-50,000 bp |---------------------------|
+ Ultra-long:     >100,000 bp      |-------------------------------------------|
```

**Why do we align?**
- To find out *where* in the genome your sequenced DNA fragment came from
- To identify differences (mutations) between your sample and the reference
- To reconstruct the complete DNA sequence from millions of short pieces

**Think of it like:**
- Having a massive jigsaw puzzle (genome)
- Getting millions of tiny pieces (reads)
- Finding where each piece fits (alignment)
- Reconstructing the full picture


#### 1. Getting the reference genome 
Before calling variants, we need to align our reads to a reference genome. We will download the reference genome for E. coli (REL606) from the NCBI website and keep it in a dedicated folder.  

Within `seq-analysis` make a dir called refgenome and navigate there: 
```diff
user1@vm-corso-colonna:~/seq-analysis$ mkdir refgenome
user1@vm-corso-colonna:~/seq-analysis$ cd refgenome/
user1@vm-corso-colonna:~/seq-analysis/refgenome$ 
```

Download, gunzip, and explore the reference genome: 
```diff 

+ COMMAND  curl - transfer a URL

user1@vm-corso-colonna:~/seq-analysis/refgenome$ curl -L -o ecoli_rel606.fasta.gz \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz

+ COMMAND gunzip 
user1@vm-corso-colonna:~/seq-analysis/refgenome$ gunzip ecoli_rel606.fasta.gz
user1@vm-corso-colonna:~/seq-analysis/refgenome$ ls 

+ check the file name, what is changed? 

user1@vm-corso-colonna:~/seq-analysis/refgenome$ cat ecoli_rel606.fasta | less 
```

#### 2. Preparing the filesystem for results 
We will make a number of directories where to put our results  
```
user1@vm-corso-colonna:~/seq-analysis$ mkdir -p results/sam results/bam results/bcf results/vcf
user1@vm-corso-colonna:~/seq-analysis$ ls 
fastqc-res  refgenome  results  trimmed  trimmed_all
user1@vm-corso-colonna:~/seq-analysis$ ls results/
bam  bcf  sam  vcf
user1@vm-corso-colonna:~/seq-analysis$ 
```

#### 3. Alignment with `bwa`
>> **[BWA](http://bio-bwa.sourceforge.net/)** (Burrows-Wheeler Aligner) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome, using the Burrows-Wheeler transform algorithm.

The reference genome must be indexed before alignment:

```diff
+ user1@vm-corso-colonna:~/seq-analysis/refgenome$ bwa index ecoli_rel606.fasta 
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

+ Why do we index? Indexing creates a searchable data structure that speeds up alignment
```

Go to `~/seq-analysis`  and align paired-end reads to the reference genome using BWA-MEM:

```diff
user1@vm-corso-colonna:~/seq-analysis$ bwa mem \
    refgenome/ecoli_rel606.fasta \
    trimmed_all/SRR2584863_1.trimmed.fastq 
    trimmed_all/SRR2584863_2.trimmed.fastq  > results/sam/SRR2584863.aligned.sam 

+ BWA MEM is optimized for reads 70bp-1Mbp
+ Outputs SAM format to stdout (we redirect to file with >)
+ Uses paired-end mode when given two input files
```

Explore the SAM file that you just produced 
```diff 
user1@vm-corso-colonna:~/seq-analysis$ cat results/sam/SRR2584863.aligned.sam | less -S 

user1@vm-corso-colonna:~/seq-analysis$ cat results/sam/SRR2584863.aligned.sam | head -3 
@SQ	SN:CP000819.1	LN:4629812
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem refgenome/ecoli_rel606.fasta trimmed_all/SRR2584863_1.trimmed.fastq trimmed_all/SRR2584863_2.trimmed.fastq
SRR2584863.1	83	CP000819.1	3339203	60	102M	=	3338771	-534	TACCGTTAACTCTCAGGATCAGGTAACCCAAAAACCCCTGCGTGACTCGGTTAAACAGGCACTGAAGAACTATTTTGCTCAACTGAATGGTCAGGATGTGAA	8(5@>3C>C@AC@@AACCCCADDBBB:DCAFFFHFJIIGEGGCCEDJJIJHHFIJGGEFDDEJIIIFGJJIIIJJJJIIJJJIJJJJIJHHHHGFFFFFCCC	NM:i:0	MD:Z:102	MC:Z:19M	AS:i:102	XS:i:0
user1@vm-corso-colonna:~/seq-analysis$ 
```

The alignment produces a file in **SAM (Sequence Alignment/Map)** format,  a tab-delimited text file containing alignment information for sequencing reads against a reference genome. The **BAM (Binary Alignment/Map)** is the compressed binary version of SAM, with reduced file size, support for indexing
and efficient random access

Key Features of `SAM/BAM` files are: 
- **Human-readable**: SAM files can be viewed and edited with text editors
- **Standardized**: Consistent format across different aligners
- **Flexible**: Supports optional fields for tool-specific information
- **Comprehensive**: Contains all alignment information in one place

The [`SAM/BAM`](https://pubmed.ncbi.nlm.nih.gov/19505943/) file contains : 
    1. **Header Section** (Optional) with metadata about the alignment (e.g. data source information, reference sequence details, algnment method/software used) 
    2. **Alignment Section** Each line represents one read alignment with standardized fields.

Each alignment line contains **11 mandatory fields** plus optional fields:
| Field | Column | Description |
|-------|--------|-------------|
| QNAME | 1 | Query/Read name |
| FLAG | 2 | Bitwise flag (alignment properties) |
| RNAME | 3 | Reference sequence name |
| POS | 4 | 1-based leftmost mapping position |
| MAPQ | 5 | Mapping quality score |
| CIGAR | 6 | CIGAR string (alignment description) |
| RNEXT | 7 | Reference name of mate/next read |
| PNEXT | 8 | Position of mate/next read |
| TLEN | 9 | Template length |
| SEQ | 10 | Read sequence |
| QUAL | 11 | ASCII of base quality scores |


#### 4. Convert SAM to BAM, sort 
>> **[SAMtools](http://www.htslib.org/)** is a suite of programs for interacting with high-throughput sequencing data in SAM/BAM format

Using the samtools program we will convert the SAM file to BAM format in three steps 
1. Go form SAM to BAM with the view command tellinmg this command that the input is in SAM format (-S) and to output BAM format (-b)
2. Sort the BAM file by coordinates (useful for next steps)
3. Index the BAM file to make it easily accessible

```diff
+ # step 1 
user1@vm-corso-colonna:~/seq-analysis$ samtools view -S -b results/sam/SRR2584863.aligned.sam > results/bam/SRR2584863.aligned.bam
user1@vm-corso-colonna:~/seq-analysis$ ls results/bam/
SRR2584863.aligned.bam

+ # step 2 
user1@vm-corso-colonna:~/seq-analysis$ samtools sort -o results/bam/SRR2584863.aligned.sorted.bam results/bam/SRR2584863.aligned.bam
user1@vm-corso-colonna:~/seq-analysis$ ls results/bam/
SRR2584863.aligned.bam  SRR2584863.aligned.sorted.bam

+ #step 3 
user1@vm-corso-colonna:~/seq-analysis$ samtools index results/bam/SRR2584863.aligned.sorted.bam
user1@vm-corso-colonna:~/seq-analysis$ ls results/bam/
SRR2584863.aligned.bam  SRR2584863.aligned.sorted.bam  SRR2584863.aligned.sorted.bam.bai
```

We can now visualize the alignment in the BAM file using samtool tview:

```diff
+ check the options for samtools tview 

user1@vm-corso-colonna:~/seq-analysis$ samtools tview 
Usage: samtools tview [options] <aln.bam> [ref.fasta]
Options:
   -d display      output as (H)tml or (C)urses or (T)ext 
   -X              include customized index file
   -p chr:pos      go directly to this position
   -s STR          display only reads from this sample or group
   -w INT          display width (with -d T only)
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
      --reference FILE
               Reference sequence FASTA FILE [null]
      --verbosity INT
               Set level of verbosity


+ visualize the alignment 
user1@vm-corso-colonna:~/seq-analysis$ samtools tview  results/bam/SRR2584863.aligned.sorted.bam refgenome/ecoli_rel606.fasta
```


#### 5. Variant calling 
>> **[FreeBayes](https://github.com/freebayes/freebayes)** is a  Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

We will use `freebayes` to produce the vcf file. First explore the options of `freebayes`,  do the variant calling, and explore the output vcf: 

```diff
user1@vm-corso-colonna:~/seq-analysis$ freebayes -h  | less 


user1@vm-corso-colonna:~/seq-analysis$ freebayes -f refgenome/ecoli_rel606.fasta results/bam/SRR2584863.aligned.sorted.bam > results/vcf/myfirst.vcf 
user1@vm-corso-colonna:~/seq-analysis$ ls results/vcf/
myfirst.vcf

+ open the file inspect it and count the number of lines 
user1@vm-corso-colonna:~/seq-analysis$ cat results/vcf/myfirst.vcf | less -S 

user1@vm-corso-colonna:~/seq-analysis$ cat results/vcf/myfirst.vcf | wc -l 
```


Many of the variants have very low quality, therefore we will filter out those variants with poor quality using bcftools: 
>> **[BCFtools](http://www.htslib.org/)** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF, including filtering, merging, comparing, and annotating variants from multiple samples.
```diff 
user1@vm-corso-colonna:~/seq-analysis$ bcftools filter -i 'QUAL>20' results/vcf/myfirst.vcf > results/vcf/myfirst.filtered.vcf 
user1@vm-corso-colonna:~/seq-analysis$ cat results/vcf/myfirst.filtered.vcf | less -S 
user1@vm-corso-colonna:~/seq-analysis$ cat results/vcf/myfirst.filtered.vcf | wc -l 
```

```diff
!EXERCISE compare the number of variants in the vcf before and after filtering
```