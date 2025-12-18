## Short read data analysis workflow



````mermaid
graph LR
    A[Sequence reads] -->|FastQC/Trimmomatic| B[FASTQ]
    B -->|BWA/Bowtie2| C[SAM/BAM]
    C -->|Clean/Sort| D[BAM]
    D -->|GATK/bcftools| E[VCF]
    
    style A fill:#FA6868,stroke:#333,stroke-width:2px,color:#607B8F
    style B fill:#5A9CB5,stroke:#333,stroke-width:2px,color:#1B3C53
    style C fill:#FACE68,stroke:#333,stroke-width:2px,color:#1B3C53
    style D fill:#FAAC68,stroke:#333,stroke-width:2px,color:#607B8F
    style E fill:#FA6868,stroke:#333,stroke-width:2px,color:#1B3C53

````



zcat SRR2584863_1.fastq.gz  | head -n 4


user1@vm-corso-colonna:~$ mkdir resres 
user1@vm-corso-colonna:~$ cd resres 
user1@vm-corso-colonna:~/resres$ fastqc -h | less 
user1@vm-corso-colonna:~/resres$ fastqc ../datadata/SRR258* -o . 


 scp  user1@212.189.205.193:/home/user1/resres/*.html  . 
