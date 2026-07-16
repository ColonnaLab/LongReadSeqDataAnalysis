# Software Used in This Course

This file lists the main software used in the lessons, with links to the official source, repository, or installation page.

## Command Line Basics

| Software | Used for | Source / installation |
|---|---|---|
| Bash / Unix shell | Command-line work, navigation, pipes, redirects, scripts | https://www.gnu.org/software/bash/ |
| GNU coreutils | Common commands such as `ls`, `pwd`, `cat`, `head`, `wc`, `mkdir` | https://www.gnu.org/software/coreutils/ |
| gzip | Compression and decompression with `gzip`, `gunzip`, `zcat` | https://www.gnu.org/software/gzip/ |
| curl | Downloading files from URLs | https://curl.se/ |
| OpenSSH | Remote login and file transfer with `ssh` and `scp` | https://www.openssh.com/ |
| less | Viewing long text files in the terminal | https://www.greenwoodsoftware.com/less/ |

## Short-Read Analysis

| Software | Used for | Source / installation |
|---|---|---|
| FastQC | Quality control reports for FASTQ files | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| Trimmomatic | Adapter and quality trimming of Illumina reads | http://www.usadellab.org/cms/?page=trimmomatic |
| Java | Running Trimmomatic `.jar` files | https://openjdk.org/ |
| BWA | Short-read alignment to a reference genome | https://github.com/lh3/bwa |
| SAMtools | SAM/BAM conversion, sorting, indexing, viewing, FASTA indexing |https://www.htslib.org/doc/samtools.html  https://github.com/samtools/samtools |
| FreeBayes | Small variant calling from aligned reads | https://github.com/freebayes/freebayes |
| BCFtools | Filtering and manipulating VCF/BCF files | https://github.com/samtools/bcftools |

## Long-Read Analysis

| Software | Used for | Source / installation |
|---|---|---|
| SRA Toolkit | Downloading SRA data with `prefetch` and `fasterq-dump` | https://github.com/ncbi/sra-tools |
| seqkit | FASTA/FASTQ summary statistics, read counts, lengths, N50 | https://github.com/shenwei356/seqkit |
| NanoPlot | Oxford Nanopore read QC plots and HTML reports | https://github.com/wdecoster/NanoPlot |
| minimap2 | Long-read alignment to a reference genome |https://lh3.github.io/minimap2/minimap2.html  https://github.com/lh3/minimap2 |
| SAMtools | BAM inspection, indexing, coverage and alignment summaries | https://github.com/samtools/samtools |
| BCFtools | VCF/BCF filtering and manipulation | https://github.com/samtools/bcftools |

## Genome Assembly

| Software | Used for | Source / installation |
|---|---|---|
| Flye | Long-read genome assembly from ONT and PacBio reads | https://github.com/mikolmogorov/Flye |
| Raven | Fast long-read genome assembly, useful for comparing assembly results | https://github.com/lbcb-sci/raven |
| miniasm | Fast overlap-based long-read assembly; usually requires polishing | https://github.com/lh3/miniasm |
| hifiasm | PacBio HiFi assembly and haplotype-resolved assembly workflows | https://github.com/chhylp123/hifiasm  https://hifiasm.readthedocs.io/en/latest/index.html|
| Canu / HiCanu | Long-read assembly, including HiFi-aware assembly modes | https://github.com/marbl/canu |
| QUAST | Assembly quality assessment and reference-based assembly comparison | https://github.com/ablab/quast |
| Bandage | Assembly graph visualization for GFA graph files | https://github.com/rrwick/Bandage |
| IGV | Visual inspection of assembly-to-reference alignments | https://igv.org/ |

## Pangenomics and Genome Graphs

| Software | Used for | Source / installation |
|---|---|---|
| PGGB | Building pangenome variation graphs from genome assemblies | https://github.com/pangenome/pggb |
| wfmash | Whole-genome pairwise alignment used by pangenome graph workflows | https://github.com/waveygang/wfmash |
| seqwish | Constructing variation graphs from all-vs-all alignments | https://github.com/ekg/seqwish |
| smoothxg | Pangenome graph normalization, smoothing, and POA consensus steps | https://github.com/pangenome/smoothxg |
| gfaffix | Graph normalization and removal of redundant graph overlaps in GFA | https://github.com/marschall-lab/GFAffix |
| ODGI | Inspecting, sorting, manipulating, and visualizing pangenome graphs | https://github.com/pangenome/odgi |
| VG | Genome graph indexing, graph read mapping with `vg giraffe`, and graph-based variant calling | https://github.com/vgteam/vg |
| PanGenie | Pangenome-based genotyping from known haplotypes and sample k-mer counts | https://github.com/eblerjana/pangenie  https://pangenie.readthedocs.io/en/latest/index.html |

## R and Reproducible Analysis

| Software | Used for | Source / installation |
|---|---|---|
| R | Statistical computing and reproducible data analysis | https://www.r-project.org/ |
| RStudio Desktop | Graphical IDE for R | https://posit.co/download/rstudio-desktop/ |

## Common Conda Installation Sources

Many bioinformatics tools above are available through Bioconda. When using conda/mamba, the usual channels are:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Useful package pages:

| Package | Bioconda page |
|---|---|
| FastQC | https://bioconda.github.io/recipes/fastqc/README.html |
| Trimmomatic | https://bioconda.github.io/recipes/trimmomatic/README.html |
| BWA | https://bioconda.github.io/recipes/bwa/README.html |
| SAMtools | https://bioconda.github.io/recipes/samtools/README.html |
| FreeBayes | https://bioconda.github.io/recipes/freebayes/README.html |
| BCFtools | https://bioconda.github.io/recipes/bcftools/README.html |
| SRA Toolkit | https://bioconda.github.io/recipes/sra-tools/README.html |
| seqkit | https://bioconda.github.io/recipes/seqkit/README.html |
| NanoPlot | https://bioconda.github.io/recipes/nanoplot/README.html |
| pycoQC | https://bioconda.github.io/recipes/pycoqc/README.html |
| minimap2 | https://bioconda.github.io/recipes/minimap2/README.html |
| Flye | https://bioconda.github.io/recipes/flye/README.html |
| Raven | https://bioconda.github.io/recipes/raven-assembler/README.html |
| miniasm | https://bioconda.github.io/recipes/miniasm/README.html |
| hifiasm | https://bioconda.github.io/recipes/hifiasm/README.html |
| Canu | https://bioconda.github.io/recipes/canu/README.html |
| QUAST | https://bioconda.github.io/recipes/quast/README.html |
| Bandage | https://bioconda.github.io/recipes/bandage/README.html |
| PGGB | https://bioconda.github.io/recipes/pggb/README.html |
| wfmash | https://bioconda.github.io/recipes/wfmash/README.html |
| seqwish | https://bioconda.github.io/recipes/seqwish/README.html |
| smoothxg | https://bioconda.github.io/recipes/smoothxg/README.html |
| gfaffix | https://bioconda.github.io/recipes/gfaffix/README.html |
| ODGI | https://bioconda.github.io/recipes/odgi/README.html |
| VG | https://bioconda.github.io/recipes/vg/README.html |
| PanGenie | https://bioconda.github.io/recipes/pangenie/README.html |
