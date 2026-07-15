# How `vg giraffe` and `vg call` Work

## Scope

This document explains two connected stages of a graph-based short-read workflow in the **vg toolkit**:

1. **`vg giraffe`** maps individual sequencing reads to a pangenome graph.
2. **`vg call`** aggregates support from mapped reads and infers genotypes for alternative paths through graph variation sites.

The explanation is based primarily on the Giraffe paper by Sirén *et al.* [1], the structural-variant genotyping paper by Hickey *et al.* [2], and the current vg documentation [3–5]. Because vg is actively developed, command-line details and internal scoring behavior can change between releases. The exact thresholds described in the 2020 `vg call` paper are therefore identified as **paper-specific algorithm details**, while the command examples follow the more recent vg documentation.

---

## 1. The overall workflow

A pangenome graph represents DNA as sequence-bearing nodes connected by edges. A path through the graph can represent a reference chromosome, an assembled haplotype, or an allele containing a SNP, insertion, deletion, inversion, or more complex combination of variants [8].

```text
FASTQ reads
    |
    v
vg giraffe
    |
    v
graph alignments (usually GAM, optionally GAF/SAM/BAM/CRAM)
    |
    v
vg pack
    |
    v
compressed node/base/edge support
    |
    v
vg call
    |
    v
sample genotypes in VCF
```

The two tools answer different questions:

- **`vg giraffe`:** Where in the pangenome did this read most likely originate?
- **`vg call`:** Given all mapped reads, which graph traversal or pair of traversals is best supported at each variation site?

A major advantage of graph mapping is that a known alternative allele is already a valid reference path. For example, a read spanning a known deletion can align directly across the deletion edge instead of being represented only as a large gap, discordant pair, or soft-clipped alignment against a single linear reference [1,2].

---

# Part I — How `vg giraffe` Works

## 2. Giraffe maps to a haplotype-informed graph

A sequence graph can contain a combinatorially large number of possible paths. Most arbitrary combinations of nearby alleles are not observed biological haplotypes. Giraffe reduces this search space by using haplotype paths stored in a **GBWT** or, in current workflows, a **GBZ** file [1,3,6].

Conceptually:

```text
Reference:           A -- B -- C -- D -- E
Deletion haplotype:  A -- B -------- D -- E
Insertion haplotype: A -- B -- X -- C -- D -- E
```

The graph encodes all branches, while the haplotype index records which combinations of branches occur together on known or sampled haplotypes. This helps Giraffe avoid treating every graph walk as equally plausible [1,6].

## 3. Main indexes used by current Giraffe workflows

A typical `vg autoindex --workflow giraffe` run produces the following files [3]:

| File | Main role |
|---|---|
| `.giraffe.gbz` | Stores the graph sequence and GBWT haplotype information |
| `.shortread.withzip.min` | Minimizer index used for seed discovery; includes zipcode annotations |
| `.dist` | Minimum-distance index for positions in the graph |
| `.shortread.zipcodes` | Stores longer zipcode information used with the minimizer/distance machinery |

The GBZ contains a GBWT haplotype index and a GBWTGraph that provides graph node sequences. The minimizer index identifies candidate locations, and the distance-related indexes help determine whether seed hits are mutually compatible in graph space [3].

## 4. Step 1: Select minimizers from the read

Giraffe uses a sparse set of representative \(k\)-mers called **minimizers** rather than querying every \(k\)-mer in a read [1].

```text
Read:        A C G T T A C G A T ...
Selected:        G T T A C
                         A C G A T
```

For each selected minimizer, the minimizer index returns occurrences along indexed graph haplotypes. Each occurrence is a **seed hit** suggesting that the read may originate near that graph position.

This is efficient because fewer strings are queried than in a full \(k\)-mer index, rare minimizers can be highly informative, and repetitive minimizers can be limited because they generate many ambiguous hits.

## 5. Step 2: Cluster seed hits using graph distance

On a linear genome, compatibility between seed hits can be evaluated using coordinate differences. A graph has branches, inversions, and multiple embedded paths, so there is no single universal coordinate axis.

Giraffe uses a distance index to group seed hits whose **minimum graph distances** are compatible with their positions in the read [1]. A promising cluster generally has seed hits in a compatible order, graph distances consistent with read-coordinate distances, and enough support to justify further alignment.

The purpose of clustering is to identify a small set of candidate graph regions before performing more expensive alignment work.

## 6. Step 3: Extend promising seeds gaplessly along haplotypes

Giraffe extends selected seed hits in both directions along haplotype-consistent graph paths, producing **gapless extensions** [1].

This is especially effective for short reads because many differences from the linear reference are already represented as graph branches:

```text
Reference path: A -- B -- C -- D
Deletion path:  A -- B ------- D

Read:           A -- B ------- D
```

From the read's perspective, there is no deletion relative to the chosen graph path. The read is simply an exact or near-exact match to an alternative haplotype.

## 7. Step 4: Perform local gapped alignment when necessary

A gapless extension may not explain an entire read because of sequencing errors, a novel SNP or indel absent from the graph, an inaccurate graph allele, or multiple seeds separated by an unresolved region.

Giraffe then performs more expensive alignment only in restricted candidate regions:

```text
seed broadly -> cluster candidates -> extend cheaply -> align expensively only where needed
```

It does not run unrestricted dynamic programming against the entire pangenome for every read [1].

## 8. Step 5: Rank mappings and calculate mapping quality

Giraffe compares candidate alignments using alignment score, seed support, consistency, and ambiguity among alternative locations. It chooses a primary alignment and assigns mapping quality according to how confidently the selected mapping is distinguished from alternatives [1].

For paired-end data, the mapper also considers read orientation, whether both ends map compatibly, estimated fragment length, and whether the pair forms a plausible fragment.

The output is commonly GAM, although current vg versions can also emit GAF or project alignments onto named reference paths for SAM/BAM/CRAM output [3].

## 9. Example Giraffe command

```bash
vg giraffe \
  -Z graph.giraffe.gbz \
  -m graph.shortread.withzip.min \
  -z graph.shortread.zipcodes \
  -d graph.dist \
  -f sample_R1.fastq.gz \
  -f sample_R2.fastq.gz \
  > sample.gam
```

In this command:

- `-Z` supplies the GBZ graph and haplotype index;
- `-m` supplies the minimizer index;
- `-z` supplies the zipcodes file;
- `-d` supplies the distance index;
- each `-f` supplies a FASTQ file [3].

---

# Part II — How `vg call` Works

## 10. What `vg call` is doing conceptually

For graph-based genotyping, `vg call` does not ask whether a read differs from one linear reference. Instead, it asks which complete paths through each graph variation site receive the strongest read support [2].

```text
mapped reads
    |
    v
coverage/support on graph bases and edges
    |
    v
variation sites defined by snarls
    |
    v
candidate traversals through each site
    |
    v
support comparison
    |
    v
genotype and VCF representation
```

The 2020 paper focused on genotyping known structural variants already represented in the graph [2].

## 11. Step 1: `vg pack` builds a compressed support index

```bash
vg pack \
  -x graph.gbz \
  -g sample.gam \
  -Q 5 \
  -o sample.pack
```

The support index records coverage over graph bases and edges. In the 2020 method, alignments with mapping quality below 5 were excluded, and the pack stored the number of qualifying reads mapped to every graph edge and every base of every node [2].

Edge support is particularly important for structural variants. For a deletion allele, reads crossing the deletion edge provide direct evidence for the breakpoint traversal, whereas node coverage inside the deleted region supports the non-deleted allele.

## 12. Step 2: Decompose the graph into snarls

A **snarl** is a graph region bounded by two node sides that encloses one or more alternative traversals. It is a graph generalization of a bubble and can represent simple, overlapping, nested, or structurally complex variation [2,7].

A simple SNP-like snarl:

```text
          /-- A --\
left ---<         >--- right
          \-- G --/
```

A deletion-like snarl:

```text
          /-- B -- C --\
left ---<              >--- right
          \------------/
```

The snarl decomposition provides bounded units within which alternative paths can be generated and compared.

## 13. Step 3: Generate candidate traversals

For each snarl, `vg call` needs a set of possible paths from the entrance boundary to the exit boundary.

When the graph was constructed from a VCF with retained allele-path information, candidate traversals correspond to combinations of VCF alleles in the snarl [2]. If a snarl contains variants \(v_1, v_2, \ldots, v_k\), and variant \(v_i\) has \(|v_i|\) alleles, the theoretical number of haplotypes is:

\[
\prod_{i=1}^{k}|v_i|
\]

This allows vg to consider overlapping variants jointly instead of always genotyping each VCF row independently.

When exact VCF allele paths are not supplied, the caller can derive candidate traversals from graph structure. In current GBZ workflows, `vg call -z` can restrict possible alleles to haplotypes represented in the GBZ, which the current documentation describes as usually faster and more accurate [5].

## 14. Step 4: Score support for each traversal

For every candidate path through a snarl, the 2020 algorithm calculated average support across the path's bases and edges using the packed coverage index [2].

```text
Reference traversal: L -> A -> B -> R
Deletion traversal:  L ------------> R
```

Evidence for the reference traversal includes coverage on `A`, `B`, and their connecting edges. Evidence for the deletion traversal includes reads supporting the direct `L -> R` edge.

The two most strongly supported paths were retained as the main diploid candidates in the 2020 algorithm [2].

## 15. Step 5: Convert traversal support into a genotype

The following rules are specific to the algorithm described in Hickey *et al.* (2020) [2]:

1. The best-supported path had to exceed a minimum support threshold, with a default of 1.
2. If it had more than \(B\) times the support of the second path, where the default \(B=6\), the site was called homozygous for the best path.
3. Otherwise, if the second path also exceeded the minimum support threshold, the site was called heterozygous for the two best paths.
4. The selected paths were translated back into VCF alleles.

```text
reference >> alternate  -> 0/0
reference ~ alternate   -> 0/1
alternate >> reference  -> 1/1
```

The paper also reported that when the number of possible haplotypes in a snarl exceeded 500,000, alleles with average support below 1 were filtered to control combinatorial growth [2].

These thresholds explain the original published algorithm but should not be assumed to describe every internal scoring detail of every current vg release.

## 16. Complete deletion example

Assume the graph contains a known deletion:

```text
Reference allele:
L -- A -- B -- C -- R

Deletion allele:
L ---------------- R
```

Reads from the non-deleted chromosome map through:

```text
L -> A -> B -> C -> R
```

Reads spanning the deletion breakpoint map through:

```text
L -> R
```

`vg pack` records support for the bases in `A`, `B`, and `C`, the reference edges, and the deletion edge `L -> R`.

| Read evidence | Expected genotype |
|---|---|
| Strong internal coverage and little `L -> R` support | `0/0` |
| Both internal and deletion-edge support | `0/1` |
| Strong `L -> R` support and little internal coverage | `1/1` |

This illustrates why Giraffe and `vg call` complement each other: Giraffe places each read on an allele-aware graph path, and `vg call` combines many mappings into a diploid genotype.

---

# Part III — Practical Usage and Important Distinctions

## 17. Current graph-native genotyping example

```bash
# 1. Map paired-end short reads
vg giraffe \
  -Z graph.giraffe.gbz \
  -m graph.shortread.withzip.min \
  -z graph.shortread.zipcodes \
  -d graph.dist \
  -f sample_R1.fastq.gz \
  -f sample_R2.fastq.gz \
  > sample.gam

# 2. Summarize read support on the graph
vg pack \
  -x graph.giraffe.gbz \
  -g sample.gam \
  -Q 5 \
  -o sample.pack

# 3. Genotype graph sites
vg call \
  graph.giraffe.gbz \
  -k sample.pack \
  -s SAMPLE \
  -z \
  > sample.vcf
```

The two uses of `-z` are unrelated:

- in `vg giraffe`, `-z` names the zipcodes index;
- in `vg call`, `-z` restricts candidate alleles to haplotypes in the GBZ [3,5].

Verify options against the installed version:

```bash
vg giraffe --help
vg pack --help
vg call --help
```

## 18. Calling graph sites versus reproducing an original VCF

### Genotype variation represented by the graph

```bash
vg call graph.gbz -k sample.pack -s SAMPLE -z
```

This emits a VCF representation of graph sites. The representation may not be identical to the VCF originally used to build the graph because conversion between VCF and graph structures can normalize, combine, or shift variants [5].

### Genotype the exact original VCF alleles

The current vg documentation states that this requires the graph to have been constructed with retained allele-path information, notably using `vg construct -a`, and that the same VCF must be supplied back to `vg call -v` [5].

A graph produced only through a normal `vg autoindex` workflow is not sufficient for this exact original-VCF mode [5].

## 19. Known-variant genotyping versus novel variant discovery

`vg call` is most straightforward when the relevant allele already exists as a graph traversal:

```text
known allele in graph
        |
reads map to it
        |
vg pack measures its support
        |
vg call genotypes it
```

If an allele is absent from the graph, reads may align with mismatches, gaps, or clipping, but the caller cannot simply select a non-existent traversal.

`vg augment` can add read-supported edits to a mutable graph before calling. However, the current vg documentation warns that augmentation-based de novo calling is not suitable for novel structural-variant discovery, and it often recommends projecting alignments to BAM and using a specialized caller for novel small variants [5].

Therefore:

- **Known SV genotyping:** graph-native Giraffe + pack + call is a natural workflow.
- **Novel SV discovery:** use a purpose-built SV discovery workflow.
- **Novel SNP/small-indel calling:** graph augmentation is possible, but specialized callers on surjected BAM may perform better [5].

## 20. Key interpretation

> **`vg giraffe` chooses the most plausible graph location and haplotype-consistent alignment for each read. `vg call` aggregates read support across graph bases and edges to choose the most plausible allele traversal or diploid pair of traversals through each variation site.**

The graph makes alternative alleles ordinary reference paths rather than merely deviations from a single privileged linear sequence. This can reduce reference bias and improve interpretation of reads supporting known structural variation [1,2,8].

---

# References

1. Sirén J, Monlong J, Chang X, Novak AM, Eizenga JM, Markello C, Sibbesen JA, Hickey G, Chang P-C, Carroll A, Gupta N, Gabriel S, Blackwell TW, Ratan A, Taylor KD, Rich SS, Rotter JI, Haussler D, Garrison E, Paten B. **Pangenomics enables genotyping of known structural variants in 5,202 diverse genomes.** *Science*. 2021;374(6574):abg8871. [https://doi.org/10.1126/science.abg8871](https://doi.org/10.1126/science.abg8871)

2. Hickey G, Heller D, Monlong J, Sibbesen JA, Sirén J, Eizenga J, Dawson ET, Garrison E, Novak AM, Paten B. **Genotyping structural variants in pangenome graphs using the vg toolkit.** *Genome Biology*. 2020;21:35. [https://doi.org/10.1186/s13059-020-1941-7](https://doi.org/10.1186/s13059-020-1941-7)

3. vgteam. **Mapping short reads with Giraffe.** vg Wiki. [https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe)

4. vgteam. **Giraffe best practices.** vg Wiki. [https://github.com/vgteam/vg/wiki/Giraffe-best-practices](https://github.com/vgteam/vg/wiki/Giraffe-best-practices)

5. vgteam. **SV Genotyping and variant calling.** vg Wiki. [https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling](https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling)

6. Sirén J, Garrison E, Novak AM, Paten B, Durbin R. **Haplotype-aware graph indexes.** *Bioinformatics*. 2020;36(2):400–407. [https://doi.org/10.1093/bioinformatics/btz575](https://doi.org/10.1093/bioinformatics/btz575)

7. Paten B, Eizenga JM, Rosen YM, Novak AM, Garrison E, Hickey G. **Superbubbles, Ultrabubbles, and Cacti.** *Journal of Computational Biology*. 2018;25(7):649–663. [https://doi.org/10.1089/cmb.2017.0251](https://doi.org/10.1089/cmb.2017.0251)

8. Garrison E, Sirén J, Novak AM, Hickey G, Eizenga JM, Dawson ET, Jones W, Garg S, Markello C, Lin MF, Paten B, Durbin R. **Variation graph toolkit improves read mapping by representing genetic variation in the reference.** *Nature Biotechnology*. 2018;36(9):875–879. [https://doi.org/10.1038/nbt.4227](https://doi.org/10.1038/nbt.4227)
