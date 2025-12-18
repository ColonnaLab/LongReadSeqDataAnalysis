markdown# Example in GitHub README.md
````mermaid
graph LR
    A[FASTQ] --> B[Quality Control]
    B --> C[Alignment]
    C --> D[Variant Calling]
````

## Types of Diagrams Mermaid Can Create

1. **Flowcharts** - Like the genomics workflow
2. **Sequence diagrams** - For showing interactions over time
3. **Gantt charts** - Project timelines
4. **Class diagrams** - Software architecture
5. **State diagrams** - System states
6. **Entity relationship diagrams** - Database schemas
7. **Git graphs** - Version control visualization

## Example for Your Genomics Workflow
````mermaid
graph LR
    A[Sequence reads Raw Data] -->|FASTQ| B[Quality ControlFastQC/Trimmomatic]
    B -->|FASTQ| C[Alignment to GenomeBWA/Bowtie2]
    C -->|SAM/BAM| D[Alignment CleanupSort/Mark Duplicates]
    D -->|BAM| E[Variant CallingGATK/bcftools]
    E -->|VCF| F[Downstream Analysis]
    
    style C fill:#ff9500,stroke:#333,stroke-width:2px,color:#fff
    style D fill:#4169e1,stroke:#333,stroke-width:2px,color:#fff
    style E fill:#ff9500,stroke:#333,stroke-width:2px,color:#fff
````

## Example for Your Genomics Workflow
````mermaid
graph LR
    A[Sequence reads Raw Data] -->|FastQC/Trimmomatic| B[FASTQ]
    B -->|BWA/Bowtie2| C[SAM/BAM]
    C -->|Clean/Sort| D[BAM]
    D -->|GATK/bcftools| E[VCF]
    
    style B fill:#5A9CB5,stroke:#333,stroke-width:2px,color:#001F3D
    style C fill:#FACE68,stroke:#333,stroke-width:2px,color:#001F3D
    style D fill:#FAAC68,stroke:#333,stroke-width:2px,color:#fff
````



zcat SRR2584863_1.fastq.gz  | head -n 4


user1@vm-corso-colonna:~$ mkdir resres 
user1@vm-corso-colonna:~$ cd resres 
user1@vm-corso-colonna:~/resres$ fastqc -h | less 
user1@vm-corso-colonna:~/resres$ fastqc ../datadata/SRR258* -o . 


 scp  user1@212.189.205.193:/home/user1/resres/*.html  . 
