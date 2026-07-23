# Tree Nine 🌲
Tree Nine is a phylogenetics workflow that forms the basis of a tuberculosis cluster tracker built with funding from the California Department of Health. Although the cluster tracker frontend is not public, this workflow can still be used to generate cluster information for any arbitrary prokaryote with a non-segmented genome. Clustering is optional; it is also possible to simply place your samples on an existing phylogenetic tree using [UShER](https://www.nature.com/articles/s41588-021-00862-7).

This repo also contains the following subworkflows:
* [Convert to Nextstrain](./convert_to_nextstrain.md) (for viewing in Auspice, non-clade sample annotations, etc)
* [matUtils annotate](./annotate.md)
* [matUtils extract](./extract.md)
* [matUtils summarize](./summarize.md)
* [Mask tree](./mask_tree.wdl)
* [Mask subtree](./mask_subtree.wdl)

## Features
* Highly scalable, even on lower-end computes
  * Preliminary development tests run directly on a seven-year-old Macbook
  * Places >11,000 new samples on a base tree of >130,000 other samples in less than two hours (on GCP, not the laptop)
* Automatic clustering (including recursive subclusters) based on genetic distance
  * Clustering can be limited to only samples specified by the user, all newly added samples, or all samples on the entire tree
  * Create per-cluster subtrees (pb/nwk)
  * Cluster IDs can be made persistent to track changes over multiple runs, but starting from scratch is also supported
* Includes a sample input tree created from SRA data if no input tree is specified 
* Trees automatically converted to UsHER (.pb), Taxonium (.jsonl.gz), Newick (.nwk), and Nextstrain (.json) formats
* Reroot the tree to a specified node
* Mimic BioNumerics's rules by "backmasking" related samples against each other to hide ambigious positions
* Summarize input, reroot, and output trees with matutils
* Filter out positions by coverage at that position and/or entire samples by overall coverage
* Annotate clades via matUtils with a specified annotation TSV

## Clustering
Clustering is an optional feature, but the main focus of the CPDH grant and ergo development. Unlike the SARS-CoV-2 cluster tracker, for tuberculosis, we determine clusters are determined based on branch length assigned by UShER and do not require locational metadata (although county-level metadata can be tracked and displayed if provided). Tree Nine's clustering script recursively searches for clusters at distances of 20, 10, and 5. Every run of Tree Nine with clustering outputs files with cluster information; if you put those files back into Tree Nine again, you can track changes to your clusters over time, maintaining persistent cluster IDs the entire time. Clusters with subclusters "know" the identity of their subclusters and vice versa, maintaining clear linkage across runs.

### Visualization of Clusters on Microreact
<img width="1536" height="876" alt="AND AGAIN" src="https://github.com/user-attachments/assets/32b6a7db-f732-41c0-a2b4-30cb08d90575" />

You can visualize your clusters on [Microreact](https://microreact.org/) with the addition of a Microreact API key, as well as share these projects with other collaborators. Tree Nine uses an undocumented beta Microreact feature to mass-share private projects with other people without requiring people to click through one email invite per project, but also supports the typical email-sharing option as a fallback. Tree Nine takes care to avoid overloading the Microreact API and has been tested uploading more than 4000 projects without issue.

If a cluster's contents change across runs, its associated Microreact project will be updated as necessary. If a cluster loses all of its samples (due to being merged, samples being reassigned, etc) its Microreact project will warn the user of this scenario, providing a list of all samples it had prior so the user can track down their missing samples. Additionally, clusters with subclusters/parent clusters link back to their parent via URL.

## Compatibility
Tree Nine takes in samples as MAPLE-formatted diff files, which you can generate with [myco](https://github.com/aofarrel/myco) or [convert from VCF](https://github.com/aofarrel/vcf_to_diff_wdl). If you are not running on tuberculosis, you will need to provide a reference genome in .fasta format, otherwise H37Rv is assumed. Any arbitrary prokaryote with a non-segmented genome should work.

For CDPH's TB cluster tracker, this system runs on [Terra](https://terra.bio/), but like any [WDL workflow](https://openwdl.org/) it can run on basically anything that supports Docker (Singularity is untested). We have verified compatibility with both [Cromwell](https://github.com/broadinstitute/cromwell) and [miniwdl](https://github.com/chanzuckerberg/miniwdl), although miniwdl requires the `--copy-input-files` runtime attribute.

Hardware values in the WDL's `runtime` sections are not minimum requirements, they are just default values we use for cloud runs. If your machine uses ARM hardware (including Apple Silicon), Docker must use compatibility layer and may have diminished performance, but it does work.
 
## Benchmarking
Formal benchmarks have not been established, but if your computer can run Docker, it can probably handle small-scale runs of the entirety of Tree Nine (placing ~10 samples on a 70 sample base tree, including clustering) in about a minute.

Placing approximately 11,000 real-world samples on a ~130,000 sample base tree takes approximately one hour, plus one hour of optional matOpimize. Clustering depends heavily on the number of samples you are considering for clustering, how large your base tree is, and how many clusters are actually found.
