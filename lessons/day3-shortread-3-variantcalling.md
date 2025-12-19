[back to course home page ](../README.md)

## Variant Calling 
This page contain a short and adapted version of the Software Carpentry lesson [Variant Calling Workflow ](https://datacarpentry.github.io/wrangling-genomics/04-variant_calling.html). We will use it as notes to key concepts we will discuss during our lesson

This workflow demonstrates how to identify genetic variants from aligned sequencing data. Each step transforms the data from raw alignments to filtered variant calls ready for downstream analysis. The pipeline uses standard tools like BWA for alignment, SAMtools for processing, and bcftools for variant calling.


#### getting the reference genome 
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

#### alignment with `bwa`

##### What Does It Mean to "Align" DNA Sequences?
**Aligning** means finding where a short piece of DNA fits best within a much longer DNA sequence - like finding where a sentence fragment belongs in a book. The computer tries millions of positions to find the best match, accounting for small differences that naturally occur between individuals.

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

+ Short Reads:    100-300 bp     |-----|
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




| variatnts | example |
|-----------|---------|
| **Mismatches** | Ref: `ATCG`<br>Read: `ATGG` |
| **Insertions** | Ref: `AT_CG`<br>Read: `ATACG` |
| **Deletions** | Ref: `ATCG`<br>Read: `AT_G` |


**[BWA](http://bio-bwa.sourceforge.net/)** (Burrows-Wheeler Aligner) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome, using the Burrows-Wheeler transform algorithm.

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
    trimmed_all/SRR2584863_2.trimmed.fastq \
    > results/sam/SRR2584863.aligned.sam 

+ BWA MEM is optimized for reads 70bp-1Mbp
+ Outputs SAM format to stdout (we redirect to file with >)
+ Uses paired-end mode when given two input files
```

##### understanding the SAM/BAM file format 
The alignment produces a file in **SAM (Sequence Alignment/Map)** format,  a tab-delimited text file containing alignment information for sequencing reads against a reference genome.

The **BAM (Binary Alignment/Map)** is the compressed binary version of SAM, with reduced file size, support for indexing
and efficient random access

Explore the SAM file that you just produced 
```diff 
user1@vm-corso-colonna:~/seq-analysis$ cat results/sam/SRR2584863.aligned.sam | less -S 

user1@vm-corso-colonna:~/seq-analysis$ cat results/sam/SRR2584863.aligned.sam | head -3 

```

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

Here is an example of a SAM entry

```
read_001  163  chr1  12345  60  100M  =  12445  200  AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG  <<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<  NM:i:1  MD:Z:99A0  AS:i:98
```

Key Features are: 
- **Human-readable**: SAM files can be viewed and edited with text editors
- **Standardized**: Consistent format across different aligners
- **Flexible**: Supports optional fields for tool-specific information
- **Comprehensive**: Contains all alignment information in one place

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



