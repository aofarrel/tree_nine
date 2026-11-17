# Extract

This workflow takes in an UShER-formatted protobuff (.pb) phylogenetic tree, and uses `matUtils extract` to either:
* extract a subtree
* collapse the tree's nodes (useful after running matOptimize)

If extracting a subtree, `samples_in_subtree` file must have the name of each sample that is to go into the subtree, one sample per line.

Note that `matUtils extract` is often used to convert a .pb into a different file format. In this repo, we instead have distinct workflows for those conversions:
* convert_to_nextstrain.wdl
* convert_to_nwk.wdl
* convert_to_taxonium.wdl