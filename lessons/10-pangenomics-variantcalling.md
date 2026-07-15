[back to course home page ](../README.md)

## Advanced 4: Pangenomics Variant Calling

This goal of this session is to use a pangenome graph for variant discovery and read-supported graph variant calling.

The practical challenge is in [the pangenomics variant-calling hackathon brief](10-pangenomics-variantcalling-hackathon.md).

Prerequisites for the practical:

- a pangenome graph built from multiple assemblies
- a list of embedded graph paths
- a candidate reference path for VCF coordinates
- optional reads for read-supported analysis

```diff
+ We will distinguish path-derived and read-supported graph variants
+ We will discuss why graph indexes are needed for read mapping
+ We will interpret graph-based VCF output relative to a chosen path
```

### **1. Variant Calling from a Pangenome Graph**

There are two related but different tasks:

- call variants directly from graph paths, using the assembled genomes embedded in the graph
- map reads to the graph and use read support to genotype or call variants

Graph/path-derived calls describe differences already represented by the assembled genome paths in the graph. Read-supported calls use sequencing reads to provide evidence for alleles in one sample.

Graph coordinates are more complex than linear reference coordinates. For VCF output, choose one graph path as the reference path and project graph variation relative to that path.

### **2. Graph Indexes and Read Mapping**

`vg giraffe` is a graph-aware read mapper. It does not call variants by itself. A typical read-supported workflow is:

```diff
+ pangenome graph -> vg indexes -> vg giraffe mapping -> vg pack -> vg call -> VCF
```

For a deeper explanation of how `vg giraffe`, `vg pack`, and `vg call` work internally, see [How `vg giraffe` and `vg call` Work](vg_giraffe_and_vg_call_explained.md).

The graph must be indexed before read mapping because graph mapping has to search across many possible paths, not only one linear reference sequence. VG Giraffe uses graph, minimizer, and distance indexes to place reads efficiently and to evaluate paths through complex regions.

### **3. Deconstruct Graph Variation from Paths**

Graph deconstruction reports variation already present in the graph paths. These calls are derived from the assembled genomes or haplotypes used to build the graph.

This is different from read-supported calling. In graph deconstruction, the graph paths themselves are the evidence. In read-supported calling, mapped reads provide evidence for alleles in a sample.

### **4. Map Reads to the Graph**

Graph-aware read mapping can reduce reference bias because reads are allowed to align to alternative alleles already represented in the graph. This is useful when the sample carries sequence that is absent from, or divergent from, a single linear reference.

Mapping quality still matters. Poor mapping can come from low-quality reads, an overly complex graph, missing sequence in the graph, repeated regions, or an inappropriate graph construction strategy.

### **5. Compute Read Support and Call Variants**

Read-supported graph variant calling summarizes how mapped reads support graph alleles. The read alignments are converted into graph coverage information, and the variant caller uses that support to produce VCF records.

VCF output still needs a coordinate system. In a graph workflow, the coordinate system is usually defined by selecting one embedded path as the reference path. Variant interpretation therefore depends on both graph construction and reference-path choice.

### **6. Key Points**

- Graph/path-derived variants come from assembled genome paths embedded in the graph
- Read-supported graph calls require read mapping, support packing, and `vg call`
- `vg giraffe` maps reads to a graph; it is not itself a variant caller
- VCF interpretation depends on the selected reference path and graph construction choices
