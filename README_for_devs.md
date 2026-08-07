# README - for devs, maintainers, and other adventurous souls

## bottlenecks
See also: https://github.com/aofarrel/tree_nine/issues/33

### distance matrix
There is a development branch that attempts to refactor the distance matrix to use a different Python package. The new distance matrix method works, but fully integrating it would require a refactor of find_clusters.py and as such is currently deprioritized. It's worth doing, and it would probably save about an hour walltime for CalTBNet in the yes-clustering case, but it is not the main bottleneck of the clustering.

### matUtils extract
Each subtree requires matUtils open and close the tree; but opening the tree takes multiple seconds. Yes, there are matUtils commands that can extract multiple subtrees at once, *but not by sample name*. I've looked into using those other extraction methods but they're not reliable for our use case; what we really need is for someone to add a function to matUtils extract that allows for extracting multiple subtrees as defined in a textfile at once. This is something that's been on the backburner for a while.

## bogus fallbacks
If you're familiar with WDL, you know WDL parsers (as a design choice of the language) do not properly understand "[iff](https://en.wikipedia.org/wiki/If_and_only_if) X happens when Y is true, and X happened, then Y is true." If you're familiar with writing complex WDLs, you additionally know that optional types (`File?` instead of `File`, etc) sometimes do not play nicely with compound types or scatter(). As a result, Tree Nine coerces some optional types into not-optionals by using select_first(), where the second value is bogus.

### example
The annotate task expects a metadata file with type File. In Tree Nine, the annotate task only runs if the metadata file, with optional type File?, is defined. 

If you're used to Normal Programming Languages you may anticipate you can input that metadata file into the annotate task like this:

```
if (defined(metadata_tsv_that_the_user_input)) {
		call annotate {
			input:
                metadata_tsv = metadata_tsv_that_the_user_input
```

Nope. Since the annotate task expects a File, not a File?, you have to use select_first(). Additionally, you can't get away with `select_first([metadata_tsv_that_the_user_input])` with most WDL parsers (some versions of miniwdl begrudgingly accept it but as usual Cromwell is fussier), so you need to define a fallback file. It doesn't really matter what the fallback file is, as long as the WDL parser thinks it's valid. 

Here, I choose the output UShER tree as the fallback. Obviously, passing in a phylogenetic tree in place of a TSV would break things -- but it won't, because barring a hypothetical bug in a WDL executor where a metadata TSV defines and then undefines itself, only the metadata TSV is ever going to be selected.

```
if (defined(metadata_tsv_that_the_user_input)) {
		call annotate {
			input:
                metadata_tsv = select_first([metadata_tsv_that_the_user_input, usher_sampled_diff.usher_tree]) # bogus fallback
```

This might beg the question as to why I don't just make the annotate task expect a File instead of a File?. The short answer is "because the one-liner select_first() hack is simpler," but also:
* Optional types get messy very quickly, especially if you are dealing with scatters.
* Having an optional File? input into a task implies that the task doesn't actually need that file to run properly -- which isn't the case here; we need a metadata TSV if we're going to annotate. This is an anti-footgun practice, but also, like, in principle, falsely tagging required files as optional makes me feel kinda gross.
