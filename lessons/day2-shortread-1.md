## Short read data analysis workflow
This workflow shows how raw DNA sequencing data is processed from initial reads to final variant calls. Each colored box represents a file format, while arrows indicate the bioinformatics tools used for transformation. The pipeline begins with FASTQ files containing raw sequencing reads and ends with VCF files containing identified genetic variants.

````mermaid
graph LR
    A[FASTQ] -->|FastQC/Trimmomatic| B[FASTQ]
    B -->|BWA/Bowtie2| C[SAM/BAM]
    C -->|Clean/Sort| D[BAM]
    D -->|GATK/bcftools| E[VCF]
    
    style A fill:#FA6868,stroke:#1B3C53,stroke-width:2px,color:#1B3C53
    style B fill:#5A9CB5,stroke:#333,stroke-width:2px,color:#1B3C53
    style C fill:#FACE68,stroke:#333,stroke-width:2px,color:#1B3C53
    style D fill:#FAAC68,stroke:#333,stroke-width:2px,color:#1B3C53
    style E fill:#FA6868,stroke:#333,stroke-width:2px,color:#1B3C53

````

```diff
+ We will go to the wrokflow while learning more bash commands: zcat,  

```
### **1. Understanding FASTQ Format**

[FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files contain raw sequencing reads with four lines per sequence:

```diff
@SRR2584863.1 HWI-ST957:244:H73TDADXX:1:1101:4712:2181/1
TTCACATCCTGACCATTCAGTTGAGCAAAATAGTTCTTCAGTGCCTGTTTAACCGAGTCACGCAGGGGTTTTTGGGTTACCTGATCCTGAGAGTTAACGGTAGAAACGGTCAGTACGTCAGAATTTACGCGTTGTTCGAACATAGTTCTG
\+
CCCFFFFFGHHHHJIJJJJIJJJIIJJJJIIIJJGFIIIJEDDFEGGJIFHHJIJJDECCGGEGIIJFHFFFACD:BBBDDACCCCAA@@CA@C>C3>@5(8&>C:9?8+89<4(:83825C(:A#########################

+ Line 1: Header with @ followed by sequence ID
+ Line 2: Raw sequence (DNA bases)
+ Line 3: Separator (+)
+ Line 4: Quality scores (ASCII encoded)
```

```diff
! EXERCISE: See where we are with $pwd  and move to  /home/user1/genomic-lesson-data/untrimmed_fastq  #pay attention to user number!! 
! see what it contains, note the gz extension of the files  

+ COMMAND: gzip, gunzip, zcat - compress or expand files
+ COMMAND: head - output the first part of files

user1@vm-corso-colonna:~/genomic-lesson-data/untrimmed_fastq$ zcat SRR2584863_1.fastq.gz  | head -n 4
@SRR2584863.1 HWI-ST957:244:H73TDADXX:1:1101:4712:2181/1
TTCACATCCTGACCATTCAGTTGAGCAAAATAGTTCTTCAGTGCCTGTTTAACCGAGTCACGCAGGGGTTTTTGGGTTACCTGATCCTGAGAGTTAACGGTAGAAACGGTCAGTACGTCAGAATTTACGCGTTGTTCGAACATAGTTCTG
\+
CCCFFFFFGHHHHJIJJJJIJJJIIJJJJIIIJJGFIIIJEDDFEGGJIFHHJIJJDECCGGEGIIJFHFFFACD:BBBDDACCCCAA@@CA@C>C3>@5(8&>C:9?8+89<4(:83825C(:A#########################

```

user1@vm-corso-colonna:~$ mkdir resres 
user1@vm-corso-colonna:~$ cd resres 
user1@vm-corso-colonna:~/resres$ fastqc -h | less 
user1@vm-corso-colonna:~/resres$ fastqc ../datadata/SRR258* -o . 


 scp  user1@212.189.205.193:/home/user1/resres/*.html  . 



### **Quality Control of Genomics Data**
This page contains a short and adapted version of the Data Carpentry lesson [Assessing Read Quality](https://datacarpentry.github.io/wrangling-genomics/02-quality-control.html). We will use it as notes to key concepts we will discuss during our lesson.

##### **Why Quality Control in Genomics?**

- **Sequencing errors** - Base calling mistakes increase towards the end of reads
- **Adapter contamination** - Sequencing adapters can appear in your reads
- **Low quality bases** - Can lead to false variant calls and misassemblies
- **Contamination** - From other organisms or samples

```diff
+ Quality control helps us identify and remove problematic data before analysis
+ Poor quality data can lead to incorrect biological conclusions
```

### **1. Understanding FASTQ Format**

FASTQ files contain raw sequencing reads with four lines per sequence:

```
@SRR2584863.1 HWI-ST957:244:H73TDADXX:1:1101:4712:2181/1
TTCACATCCTGACCATTCAGTTGAGCAAAATAGTTCTTCAGTGCCTGTTTAACCGAGTCACGCAGGGGTTTTTGGGTTACCTGATCCTGAGAGTTAACGGTAGAAACGGTCAGTACGTCAGAATTTACGCGTTGTTCGAACATAGTTCTG
+
CCCFFFFFGHHHHJIJJJJIJJJIIJJJJIIIJJGFIIIJEDDFEGGJIFHHJIJJDECCGGEGIIJFHFFFACD:BBBDDACCCCAA@@CA@C>C3>@5(8&>C:9?8+89<4(:83825C(:A#########################
```

```diff
+ Line 1: Header with @ followed by sequence ID
+ Line 2: Raw sequence (DNA bases)
+ Line 3: Separator (+)
+ Line 4: Quality scores (ASCII encoded)
```

##### **Quality Score Encoding**

- Phred quality scores: Q = -10 × log10(P), where P is error probability
- Q30 = 99.9% accuracy (1 error in 1000 bases)
- Q20 = 99% accuracy (1 error in 100 bases)

### **2. Running FastQC**

```diff
+ COMMAND fastqc performs quality control checks on raw sequence data
```

```bash
user1@vm-corso-colonna:~/dc_workshop/data/untrimmed_fastq$ fastqc SRR2584863_1.fastq.gz

Started analysis of SRR2584863_1.fastq.gz
Approx 5% complete for SRR2584863_1.fastq.gz
Approx 10% complete for SRR2584863_1.fastq.gz
...
Analysis complete for SRR2584863_1.fastq.gz
```

##### **Running FastQC on Multiple Files**

```diff
+ Use wildcards to process multiple files at once
```

```bash
user1@vm-corso-colonna:~/dc_workshop/data/untrimmed_fastq$ fastqc *.fastq*
```

##### **FastQC Output Files**

```bash
user1@vm-corso-colonna:~/dc_workshop/data/untrimmed_fastq$ ls
SRR2584863_1_fastqc.html  SRR2584863_1_fastqc.zip  SRR2584863_1.fastq.gz
SRR2584863_2_fastqc.html  SRR2584863_2_fastqc.zip  SRR2584863_2.fastq.gz
```

### **3. Interpreting FastQC Results**

##### **Key Quality Metrics**

```diff
+ Per base sequence quality - Shows quality scores across all bases
+ Per sequence quality scores - Distribution of mean quality per read  
+ Per base sequence content - Nucleotide proportions at each position
+ Sequence length distribution - Read length consistency
+ Adapter content - Presence of sequencing adapters
```

##### **Understanding Quality Plots**

- **Green** = Good quality
- **Orange** = Warning (acceptable)
- **Red** = Poor quality (needs attention)

```diff
! EXERCISE: Open the FastQC HTML report and identify which metrics pass/warn/fail
```

### **4. Quality Trimming with Trimmomatic**

```diff
+ COMMAND trimmomatic removes adapters and trims low quality bases
```

##### **Basic Trimmomatic Command Structure**

```bash
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
```

##### **Common Trimming Options**

```diff
+ ILLUMINACLIP - Remove adapters
+ SLIDINGWINDOW - Scan with sliding window, trim when quality drops
+ LEADING - Trim low quality bases from start
+ TRAILING - Trim low quality bases from end
+ MINLEN - Drop reads shorter than specified length
```

##### **Example Trimming Command**

```bash
user1@vm-corso-colonna:~/dc_workshop/data/untrimmed_fastq$ trimmomatic PE \
  SRR2584863_1.fastq.gz SRR2584863_2.fastq.gz \
  SRR2584863_1.trim.fastq.gz SRR2584863_1un.trim.fastq.gz \
  SRR2584863_2.trim.fastq.gz SRR2584863_2un.trim.fastq.gz \
  SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

Input Read Pairs: 1107090
Both Surviving: 885220 (79.96%)
Forward Only Surviving: 216472 (19.55%)
Reverse Only Surviving: 2850 (0.26%)
Dropped: 2548 (0.23%)
```

### **5. Automating Quality Control**

##### **Creating a Trimming Script**

```bash
user1@vm-corso-colonna:~/dc_workshop/scripts$ nano trimming_script.sh
```

```bash
#!/bin/bash

# Script for automating quality trimming

for infile in *_1.fastq.gz
do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE ${infile} ${base}_2.fastq.gz \
                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
done
```

```diff
+ Make script executable and run
```

```bash
user1@vm-corso-colonna:~/dc_workshop/scripts$ chmod +x trimming_script.sh
user1@vm-corso-colonna:~/dc_workshop/scripts$ ./trimming_script.sh
```

### **6. Post-Trimming Quality Check**

```diff
! EXERCISE: Run FastQC on trimmed files and compare with original reports
```

```bash
user1@vm-corso-colonna:~/dc_workshop/data/trimmed_fastq$ fastqc *.trim.fastq.gz
user1@vm-corso-colonna:~/dc_workshop/data/trimmed_fastq$ ls *.html
SRR2584863_1.trim_fastqc.html  SRR2584863_2.trim_fastqc.html
```

##### **Key Improvements After Trimming**

- Higher average quality scores
- More uniform sequence length
- Reduced adapter content
- Better per-base sequence content

### **7. MultiQC for Summary Reports**

```diff
+ COMMAND multiqc aggregates results from multiple tools into a single report
```

```bash
user1@vm-corso-colonna:~/dc_workshop/results$ multiqc .

[INFO   ]         multiqc : This is MultiQC v1.9
[INFO   ]         multiqc : Searching : .
[INFO   ]          fastqc : Found 8 reports
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
```

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