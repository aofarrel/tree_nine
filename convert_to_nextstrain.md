# convert to nextstrain

Converts an UShER-formatted pb into an Auspice-compatiable tree or subtrees in Nextstrain format.

Metadata files are optional, but if included, need the following format (paraphrased from `matUtils extract --help`)
* TSVs or CSVs
* Sample identifiers in the first column and an arbitrary number of metadata values in separate columns
* Each file needs a header line