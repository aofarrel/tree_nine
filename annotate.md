# Annotate

Annotate a .pb with clade-specific metadata via `matutils annotate`. This should not be confused with Nextstrain-style annotation. Please see [the UsHER wiki's page on annotation](https://usher-wiki.readthedocs.io/en/latest/matUtils.html#annotate) for more information.

If is_mutation_paths (default: false), it is assumed your TSV is clade paths, and matUtils will be run with -P
'''
--clade-paths (-P): Path to a tsv file mapping clades to mutation paths which must exist in the tree.  Format is the same as the first and third columns of the output of matUtils extract --clade-paths.
'''

Otherwise, it is assumed your TSV is a simple clade-tab-samplename format, and matUtils will be run with -c
'''
--clade-names (-c): Path to a file containing clade asssignments of samples. An algorithm automatically locates and annotates clade root nodes.
'''