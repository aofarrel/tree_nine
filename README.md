# Tree Nine
Put diff files on an existing phylogenetic tree using [UShER](https://www.nature.com/articles/s41588-021-00862-7)'s `usher sampled` task with a bit of help from [SRANWRP](https://www.github.com/aofarrel/SRANWRP), followed by conversion of that tree to Taxonium, Newick, and Nextstrain formats. Samples' SNP distance is calculated and output as a distance matrix, and samples will be placed into clusters based on the distance.

Verified on Terra-Cromwell and miniwdl. Make sure to add `--copy-input-files` for miniwdl. Default inputs assume you're working with _Mycobacterium tuberculosis_, be sure to change them if you aren't working with that bacterium.

This repo also contains the following subworkflows:
* [Annotate](./annotate.md)
* [Convert to Nextstrain](./convert_to_nextstrain.md) (for viewing in Auspice, non-clade sample annotations, etc)
* [Extract](./extract.md)
* [Mask tree](./mask_tree.wdl)
* [Mask subtree](./mask_subtree.wdl)
* [Summarize](./summarize.md)

## features
* Highly scalable, even on lower-end computes
* Can input a single pre-combined diff file 
* Includes a sample input tree created from SRA data if no input tree is specified 
* Trees automatically converted to UsHER (.pb), Taxonium (.jsonl.gz), Newick (.nwk), and Nextstrain (.json) formats
* Automatic clustering based on configurable genetic distance
  * Nextstrain tree(s) will be annotated by cluster
  * Clustering can be limited to only samples specified by the user, all newly added samples, or all samples
  * Clustering is also performed after backmasking
  * (optional) Create per-cluster Nextstrain subtrees
* (optional) Reroot the tree to a specified node
* (optional) Backmask newly-added samples against each other to hide positions where any newly-added sample lacks data, then create a new set of trees based on the backmasked diff files
  * Designed for highly clonal samples which have a plausible direct epidemiological relationship 
  * Backmasking can only be performed on samples which have a sample-level diff files
* (optional) Summarize input, reroot, and output trees with matutils
* (optional) Filter out positions by coverage at that position and/or entire samples by overall coverage
* (optional) Specify your own reference genome if you don't want to work with H37Rv
* (optional) Annotate clades via matutils with a specified annotation TSV
 
## benchmarking
Formal benchmarks have not been established, but a full run of placing 60 new TB samples on an existing 7000+ TB sample tree, conversion to taxonium and newick formats, distance matrixing, clustering finding, and creating cluster-specific Nextstrain trees executes in about five minutes on a 2019 Macbook Pro.

Backmasking is the least scalable part of the pipeline. The comparison itself theoretically scales <i>n<sup>2</sup></i> and once the comparison is completed, <i>n</i> backmasked disk files must be written to the disk. We have observed that memory problems tend to arise during the file-writing part when <i>n≈55</i> on a local machine. Runtime attributes are adjustable as task-level variables to aid with scaling on cloud backends, although we have seen the default handle 60 samples at a time without much issue.