# Mask-Tree / Mask-Subtree

These workflows use `matUtils mask` to mask a given set of positions in an UShER-formatted (.pb) phylogenetic tree.

The basic masking workflow takes in the `pb_tree` and a file listing which samples to mask. This mask file needs to be formatted like this:
```
N10000002N
N10000003N
N10000009N
N12345000N
```
Where the first `N` represents the reference allele, the number in the middle represents the position, and the final `N` represents the allele to mask. You do not need to enter actual alleles -- if this was input as-is, any SNP at 10000002 will be masked, as well as 10000003, 10000009, and 12345000.

Mask-Subtree takes in an additional file representing which samples should go into the subtree. These samples should be in the file, one sample name per line.