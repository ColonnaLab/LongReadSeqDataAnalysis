


````mermaid
graph TD
    A[Sequence reads] -->|FASTQ| B[Quality control]
    B -->|FASTQ| C[Alignment to Genome]
    C -->|SAM/BAM| D[Alignment Cleanup]
    D -->|BAM| E[Variant Calling]
    E -->|VCF| F[...]
    
    style C fill:#ff9500,stroke:#333,stroke-width:2px
    style D fill:#4169e1,stroke:#333,stroke-width:2px
    style E fill:#ff9500,stroke:#333,stroke-width:2px
```

https://en.wikipedia.org/wiki/FASTQ_format

zcat SRR2584863_1.fastq.gz  | head -n 4


user1@vm-corso-colonna:~$ mkdir resres 
user1@vm-corso-colonna:~$ cd resres 
user1@vm-corso-colonna:~/resres$ fastqc -h | less 
user1@vm-corso-colonna:~/resres$ fastqc ../datadata/SRR258* -o . 


 scp  user1@212.189.205.193:/home/user1/resres/*.html  . 
