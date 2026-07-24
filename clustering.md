# Clustering
Clustering is done recursively at genetic distances 20, 10, and 5. These are sometimes informally called "20 SNPs" etc but strictly speaking they can result from indels or coverage differences. For the time being, these values are hardcoded, but the way clustering works means you can easily filter/ignore values at distances not valid for your analysis -- for example, some users of the TB Cluster Tracker ignore clusters at 20 as that's a bit excessive for tuberculosis.

Clustering is specifically defined by "every sample in a cluster is with X distance of at least one other sample in a cluster." This means that any two samples might be >X apart, as long as there exists another sample "chaining" them together. For example, when clustering at 20, if A:B = 15, B:C = 6, and A:C = 21, all three will form a cluster at 20 due to B "chaining" A and C together. Additionally, B and C will form a subcluster at 10, but that subcluster will not contain A.

A sample can only be in one cluster *at a given distance* at a time. In the persistent case (see below), if a new sample would bridge the gap between two existing clusters, those two clusters will merge into one.

## Persistent Clusters
Tree Nine is used by CDPH to track TB clusters over time by taking in inputs from previous Tree Nine runs. The recommended way of doing this is, is to have every run of Tree Nine take in:
* The exact same base tree every time (that is to say: do not use the tree output by previous Tree Nine runs)
* Last run's combined diff file
* Last run's sample names file
* Last run's persistentIDSYYYY-MM-DD.tsv file
* Last run's persistentMETAYYYY-MM-DD.tsv file
* Last run's all_cluster_informationYYYY-MM-DD.json file
* Your latest "batch" of diff files will be input as `diffs` with WDL type `Array[File]`

## Backmasking
UShER's assignment of branch length when it comes to ambigious (low-coverage) positions can vary slightly (typically only +/- 1) across runs. Other phylogenetics programs such as the (now-defunct) BioNumerics software instead do not count any ambigious positions within closely related samples. Tree Nine attempts to take the best of both worlds by identifying clusters via UShER's assignments, and within those clusters, calculates a distance matrix and subtree using both UShER's determinations and the mask-everything-ambigious rules, the latter of which we call "backmasking" (since you're going back and masking).

### Meddling with persistent cluster files
Every cluster referenced in the persistentIDs file should also be referenced in the persistentMETA file.

#### Custom cluster names
Although find_clusters.py assigns clusters a zfilled six-digit UUID, nothing is stopping you from providing a persistentIDS/persistentMETA file with arbitrary strings as cluster names -- in fact, this is what our test set does to demonstrate what's being tested. However, the perl script in process_clusters.py uses grep in such a way that persistent cluster IDs cannot be substrings of other persistent IDs. For example, if a persistent ID is called "ALPHA_BETA" and another persistent ID is called "ALPHA", that will cause problems. However, "ALPHA_BETA" and "AAAAA" would be fine, even though both share the letter A. This **only** happens if you're naming clusters outside of Tree Nine; cluster names assigned by find_clusters.py cannot generate substrings as they are [zfilled](https://www.w3schools.com/Python/ref_string_zfill.asp).

### Meddling with all_clusters_informationYYYY-MM-DD.json
Don't.

## Bottlenecks
The two largest bottlenecks in clustering are:
* The creation of the initial distance matrix between all sample considered for clustering
  * [PhyloDM](https://github.com/aaronmussig/PhyloDM) is being considered as a replacement
* Generating cluster subtrees; due to how matUtils works, this requires opening and closing the base tree once per subtree
  * The base tree used for CDPH's tuberculosis clustering takes approximately 2 seconds to open, which is hella quick all things considered, but when you're generating 8000 subtrees adds up quickly
 
Any future attempts to make clustering more efficient should focus on those two problems.
