# Tree Nine
The Tree Nine system is the basis of a tuberculosis cluster tracker built with funding from the California Department of Health, but it can also be run as a standalone workflow to simply place samples on an existing phylogenetic tree using [UShER](https://www.nature.com/articles/s41588-021-00862-7), followed by conversion of the resulting tree to newick (nwk), Taxonium, and Nextstrain (Auspice) formats. For CDPH's TB cluster tracker, this system runs on Terra, but like any WDL workflow it can run on basically anything that supports Docker (Singularity is untested). We have verified compataibility with both Cromwell and miniwdl.

If you enable clustering, samples' SNP distance from each other is calculated from branch length, allowing them to be placed into clusters. Clustering can be done on a subset of samples or across the entire tree. Every run of Tree Nine with clustering outputs files with cluster information; if you put those files back into Tree Nine again, you can track changes to your clusters over time, maintaining persistent cluster IDs the entire time. Every cluster can optionally generate a Microreact template JSON and communicate with the Microreact API to automatically create projects containing the cluster's subtree, distance matrix, links to parent/subclusters, and a sample-level metadata table.

Tree Nine takes in samples as MAPLE-formatted diff files, which you can generate with [myco](https://github.com/aofarrel/myco) or [convert from VCF](https://github.com/aofarrel/vcf_to_diff_wdl).

This repo also contains the following subworkflows:
* [Annotate](./annotate.md)
* [Convert to Nextstrain](./convert_to_nextstrain.md) (for viewing in Auspice, non-clade sample annotations, etc)
* [Extract](./extract.md)
* [Mask tree](./mask_tree.wdl)
* [Mask subtree](./mask_subtree.wdl)
* [Summarize](./summarize.md)
* Python scripts relating to clustering (run in the context of clustering)

## features
* Highly scalable, even on lower-end computes
  * Preliminary development tests run directly on a seven-year-old Macbook in about three minutes (don't worry, we test it on GCP too)
* Here's One I Made Earlier: You can input a single pre-combined diff file, or even build a continuous combined diff file over multiple runs, to keep analyzing the same subset of samples again and again to detect possible transmission clusters
* Includes a sample input tree created from SRA data if no input tree is specified 
* Trees automatically converted to UsHER (.pb), Taxonium (.jsonl.gz), Newick (.nwk), and Nextstrain (.json) formats
* Automatic clustering based on configurable genetic distance
  * Nextstrain tree(s) will be annotated by cluster
  * Clustering can be limited to only samples specified by the user, all newly added samples, or all samples on the entire tree
  * Create per-cluster subtrees (pb/nwk/Nextstrain)
  * Cluster IDs can be made persistent across runs, but starting from scratch is also supported
* Reroot the tree to a specified node
* Mimic BioNumerics's rules by "backmasking" related samples against each other to hide ambigious positions
  * Designed for highly clonal samples which have a plausible direct epidemiological relationship
  * Backmasking can only be performed on samples which have a sample-level diff files
* Summarize input, reroot, and output trees with matutils
* Filter out positions by coverage at that position and/or entire samples by overall coverage
* Specify your own reference genome if you don't want to work with H37Rv
* Annotate clades via matutils with a specified annotation TSV

## compatibility
Hardware values in the WDL's `runtime` sections are not minimum requirements, they are just default values we use for cloud runs. If your machine uses ARM hardware (including Apple Silicon), Docker must use compatiability layer and may have diminished performance, but it does work. Singularity is reported to work with [myco](https://github.com/aofarrel/myco) with some adjustments, but Singularity not been tested with Tree Nine. If running this workflow with miniwdl, include the `--copy-input-files` runtime attribute.
 
## benchmarking
Formal benchmarks have not been established, but if your computer can run Docker, it can probably handle small-scale runs of the entirety of Tree Nine (placing ~10 samples on a 70 sample base tree, including clustering) in about a minute.

Placing approximately 11,000 real-world samples on a ~130,000 sample base tree takes approximately one hour, plus one hour of optional matOpimize. Clustering depends heavily on the number of samples you are considering for clustering and how many clusters are actually found. Due to a matUtils limitation, we currently must open and close the tree file multiple times per cluster (which only takes ~2 seconds on >130,000 sample trees, but that adds up quickly if you have >3500 clusters).
