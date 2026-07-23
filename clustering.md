# Clustering
Clustering is done recursively at 20, 10, and 5. For the time being, these values are hardcoded.

## Backmasking
UShER's assignment of branch length when it comes to ambigious (low-coverage) positions can vary slightly (typically only +/- 1) across runs. Other phylogenetics programs such as the (now-defunct) BioNumerics software instead do not count any ambigious positions within closely related samples. Tree Nine attempts to take the best of both worlds by identifying clusters via UShER's assignments, and within those clusters, calculates a distance matrix and subtree using both UShER's determinations and the mask-everything-ambigious rules, the latter of which we call "backmasking" (since you're going back and masking).

## Meddling with persistent cluster files
Every cluster referenced in the persistentIDs file should also be referenced in the persistentMETA file.

### Custom cluster names
Although find_clusters.py assigns clusters a zfilled six-digit UUID, nothing is stopping you from providing a persistentIDS/persistentMETA file with arbitrary strings as cluster names -- in fact, this is what our test set does to demonstrate what's being tested. However, the perl script in process_clusters.py uses grep in such a way that persistent cluster IDs cannot be substrings of other persistent IDs. For example, if a persistent ID is called "ALPHA_BETA" and another persistent ID is called "ALPHA", that will cause problems. However, "ALPHA_BETA" and "AAAAA" would be fine, even though both share the letter A. This **only** happens if you're naming clusters outside of Tree Nine; cluster names assigned by find_clusters.py cannot generate substrings as they are [zfilled](https://www.w3schools.com/Python/ref_string_zfill.asp).

## Bottlenecks
The two largest bottlenecks in clustering are:
* The creation of the initial distance matrix between all sample considered for clustering
  * [PhyloDM](https://github.com/aaronmussig/PhyloDM) is being considered as a replacement
* Generating cluster subtrees; due to how matUtils works, this requires opening and closing the base tree once per subtree
  * The base tree used for CDPH's tuberculosis clustering takes approximately 2 seconds to open, which is hella quick all things considered, but when you're generating 8000 subtrees adds up quickly
 
Any future attempts to make clustering more efficient should focus on those two problems.
