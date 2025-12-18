

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