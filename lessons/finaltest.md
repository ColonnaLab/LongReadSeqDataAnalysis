I'll create a quiz based on these shell and bioinformatics lessons. Here's a 15-question assessment that covers the key concepts:

# Shell & Bioinformatics Quiz

**Name: ________________Date: ________________**
*Instructions: Answer all questions. Write commands exactly as you would type them.*
---

## Part 1: Basic Shell Commands (Questions 1-6)

**1.** What command would you use to find out which directory you are currently in?

Answer: ______________________

**2.** You are in `/home/user1/` and want to navigate to `/home/user1/seq-analysis/results/`. Write the command:

Answer: ______________________

**3.** What is the difference between these two commands? (2 points)
- `ls -s`
- `ls -S`

Answer: ______________________

**4.** Match the following commands with their functions:
- `mkdir` ___    a) Copy files
- `cp` ___       b) List directory contents  
- `mv` ___       c) Create directories
- `ls` ___       d) Move/rename files

**5.** What does the `>` symbol do in this command?
```bash
wc -l *.pdb > lengths.txt
```

Answer: ______________________

**6.** What is the purpose of the pipe symbol `|` in shell commands? Give an example.

Answer: ______________________

---

## Part 2: File Formats & Quality Control (Questions 7-10)

**7.** In a FASTQ file, what does each of the 4 lines represent for a single sequence? (4 points)
- Line 1: ______________________
- Line 2: ______________________
- Line 3: ______________________
- Line 4: ______________________

**8.** What does a Phred quality score of Q30 mean?
- [ ] 90% accuracy
- [ ] 99% accuracy  
- [ ] 99.9% accuracy
- [ ] 99.99% accuracy

**9.** To view the contents of a compressed file `sample.fastq.gz` without uncompressing it, which command would you use?

Answer: ______________________

**10.** In FastQC results, what do these colors indicate?
- Green: ______________________
- Orange: ______________________
- Red: ______________________

---

## Part 3: Sequence Analysis Workflow (Questions 11-15)

**11.** Put these file formats in the correct order for a variant calling workflow:
___ VCF    ___ FASTQ    ___ SAM    ___ BAM

**12.** Why must we index the reference genome before alignment?

Answer: ______________________

**13.** What is the main difference between SAM and BAM file formats? (2 points)

Answer: ______________________

**14.** Complete this command to align paired-end reads using BWA:
```bash
bwa mem refgenome/ecoli_rel606.fasta \
    _________________ \
    _________________ \
    > results/sam/aligned.sam
```

**15.** After variant calling with freebayes, you have 500 variants. You filter with:
```bash
bcftools filter -i 'QUAL>20' input.vcf > output.vcf
```
What does this command do and why is it important?

Answer: ______________________

---

**Bonus Question (+2 points):** 
What does SSH stand for and why do bioinformaticians need to use it?

Answer: ______________________

---

**Total Score: ___/20**

---

## Answer Key (For Instructor)

1. `pwd`
2. `cd seq-analysis/results/` or `cd ~/seq-analysis/results/`
3. `-s` shows file sizes, `-S` sorts by size
4. c, a, d, b
5. Redirects output to a file instead of screen
6. Passes output from one command as input to another; Example: `sort file.txt | head -n 5`
7. Header/ID, DNA sequence, Separator (+), Quality scores
8. 99.9% accuracy
9. `zcat sample.fastq.gz`
10. Green = Good quality, Orange = Warning, Red = Poor quality
11. 1-FASTQ, 2-SAM, 3-BAM, 4-VCF
12. Creates searchable data structure for faster alignment
13. SAM is text/human-readable, BAM is binary/compressed
14. `trimmed_all/SRR2584863_1.trimmed.fastq` and `trimmed_all/SRR2584863_2.trimmed.fastq`
15. Filters variants keeping only those with quality score >20; removes low-quality false positive variants

Bonus: Secure Shell; needed to connect to remote servers/computing clusters for analysis