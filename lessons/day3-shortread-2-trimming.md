[back to course home page ](../README.md)

## Quality Trimming with Trimmomatic
This page contain a short and adapted version of the Software Carpentry lesson [Assessing Read Quality](https://datacarpentry.github.io/wrangling-genomics/02-quality-control.html). We will use it as notes to key concepts we will discuss during our lesson




**[`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic)** is a flexible read trimming tool for Illumina NGS data that performs quality filtering and adapter removal to clean raw sequencing reads before downstream analysis.


##### **Basic Trimmomatic Command Structure**

```diff
trimmomatic PE \
  -threads 4 \
  -phred33 \
  input_forward.fq.gz \
  input_reverse.fq.gz \
  output_forward_paired.fq.gz \
  output_forward_unpaired.fq.gz \
  output_reverse_paired.fq.gz \
  output_reverse_unpaired.fq.gz \
  OPTION:VALUE...

# Common Trimming Options
+ ILLUMINACLIP - Remove adapters
+ SLIDINGWINDOW - Scan with sliding window, trim when quality drops
+ LEADING - Trim low quality bases from start
+ TRAILING - Trim low quality bases from end
+ MINLEN - Drop reads shorter than specified length
```

We will apply trimmomatic to our data. First create a folder in seq-analysis called `trimmed`: 

```bash

user1@vm-corso-colonna:~$  #start from your home 

user1@vm-corso-colonna:~$ cd seq-analysis/  #navigate to seq-analysis 
user1@vm-corso-colonna:~/seq-analysis$ 

user1@vm-corso-colonna:~/seq-analysis$ mkdir trimmed  #make a folder named trimmed 

user1@vm-corso-colonna:~/seq-analysis$ cd trimmed # move to trimmed 

user1@vm-corso-colonna:~/seq-analysis$ java -jar /data/Trimmomatic-0.39/trimmomatic-0.39.jar \
  PE -threads 4 \
  ../../genomic-lesson-data/untrimmed_fastq/SRR2584863_1.fastq.gz \
  ../../genomic-lesson-data/untrimmed_fastq/SRR2584863_2.fastq.gz \
  SRR2584863_1.trimmed.fastq SRR2584863_1.un.trimmed.fastq \
  SRR2584863_2.trimmed.fastq SRR2584863_2.un.trimmed.fastq \
  ILLUMINACLIP:/data/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:4:20
```

##### Understanding trimmomatic output 


```bash 
user1@vm-corso-colonna:~/seq-analysis/trimmed$ ls -l 
total 962376
-rw-rw-r-- 1 user1 user1 513256525 Dec 18 22:19 SRR2584863_1.trimmed.fastq
-rw-rw-r-- 1 user1 user1  20448026 Dec 18 22:19 SRR2584863_1.un.trimmed.fastq
-rw-rw-r-- 1 user1 user1 451612545 Dec 18 22:19 SRR2584863_2.trimmed.fastq
-rw-rw-r-- 1 user1 user1    139890 Dec 18 22:19 SRR2584863_2.un.trimmed.fastq

```
Count lines and compare with Untrimmed 

We have already trimmed all in trimmed_all 

```diff
! EXERCISE: Run FastQC on trimmed files and compare with original reports
```

```bash

```

##### **Key Improvements After Trimming**

- Higher average quality scores
- More uniform sequence length
- Reduced adapter content
- Better per-base sequence content


##### **Quality Control Best Practices**

1. **Always check raw data quality first**
2. **Document trimming parameters used**
3. **Compare before/after trimming metrics**
4. **Keep untrimmed files as backup**
5. **Check for consistent quality across samples**

```diff
+ Poor quality data → Poor results
+ Time spent on QC → Time saved debugging later
```