# Tree Nine
Master repo for the following:
* The eponymous Tree Nine workflow: A WDL workflow for placing tuberculosis samples on an existing phylogenetic tree with UShER, converting that tree into multiple output formats, then clustering said samples based on genetic distance **<--- if you are a public health department, this is probably what you're looking for**
* WDLizations of several common matUtils commands
  * matUtils annotate
  * matUtils extract (for subtrees and for file conversion)
  * matUtils mask
  * matUtils summarize
* WDLization of [UShER](https://www.nature.com/articles/s41588-021-00862-7)

For organizational sake, WDLizations of matUtils and/or UShER may be split off into different repo in the future.

## Getting Tree Nine running
Tree Nine is verifed on Terra-Cromwell and miniwdl. Make sure to add `--copy-input-files` for miniwdl. **The clustering portion of the workflow may not work on ARM machines (including Apple Silicon) due to a known Rosetta incompatibility (everything else should be fine though).**

## Tree Nine features
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
* (optional) "Backmask" newly-added samples against each other to hide positions where any newly-added sample lacks data, then create a new set of trees based on the backmasked diff files
  * Designed for highly clonal samples which have a plausible direct epidemiological relationship 
  * Backmasking can only be performed on samples which have a sample-level diff files
* (optional) Summarize input, reroot, and output trees with matutils
* (optional) Filter out positions by coverage at that position and/or entire samples by overall coverage
* (optional) Specify your own reference genome if you don't want to work with H37Rv
* (optional) Annotate clades via matutils with a specified annotation TSV
* (optional) Upload to clusters to Microreact for visualization
 
## benchmarking
Formal benchmarks have not been established, but a full run of placing 60 new TB samples on an existing 7000+ TB sample tree, conversion to taxonium and newick formats, distance matrixing, clustering, and creating cluster-specific Nextstrain trees executes in about five minutes on a 2019 Macbook Pro.