import os
import argparse
import logging
from datetime import date
import time
import numpy as np
import ete3
import polars as pl # this is overkill and takes forever to import; too bad!
pl.Config.set_tbl_rows(200)


parser = argparse.ArgumentParser(description="Crunch data, extract trees, upload to MR, etc")
parser.add_argument('-cs', '--currentsamples', type=str, help='TSV: current sample information')
#parser.add_argument('-sm', '--samplemeta', type=str, help='TSV: sample metadata pulled from terra (including myco outs), one line per sample')
#parser.add_argument('-pcm', '--persistentclustermeta', type=str, help='TSV: persistent cluster metadata from last full run of TB-D')
#parser.add_argument('-ccm', '--currentclustermeta', type=str, help='TSV: current cluster metadata (as identified by find_clusters.py like five minutes ago)')
parser.add_argument('-pid', '--persistentids', type=str, help='TSV: persistent IDs from last full run of TB-D')
#parser.add_argument('-cid', '--currentids', type=str, help='TSV: current cluster IDs (as identified by find_clusters.py like five minutes ago)')
#parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')

args = parser.parse_args()
current_samples = pl.read_csv(args.currentsamples, separator="\t")
persistent_samples_20 = pl.read_csv(args.persistent_samples_20, )
#current_clusters = pl.read_csv(args.currentclustermeta, sep="\t")
#persistent_clusters = pl.read_csv(args.persistentclustermeta, sep="\t")


# ensure each sample in current-clusters has, at most, one 20 SNP, one 10 SNP, and one 05 SNP
check_clusters_valid = current_samples.group_by("sample_id", maintain_order=True).agg(pl.col("cluster_distance"))
for row in grouped_by_sample_id.iter_rows(named=True):
    if len(row["cluster_distance"]) == 1:
        if row["cluster_distance"] == [-1]: # unclustered -- but not currently in find_clusters.py
            pass
        else:
            assert row["cluster_distance"] == [20], f"{row['sample_id']} has one cluster but it's not 20 SNP: {row['cluster_distance']}"
    elif len(row["cluster_distance"]) == 2:
        assert row["cluster_distance"] == [20,10], f"{row['sample_id']} has two clusters but it's not 20-10: {row['cluster_distance']}"
    elif len(row["cluster_distance"]) == 3:
        assert row["cluster_distance"] == [20,10,5], f"{row['sample_id']} has three clusters but it's not 20-10-5: {row['cluster_distance']}"
    else:
        print(f"{row['sample_id']} has invalid clusters: {row['cluster_distance']}")
        raise ValueError

# output current cluster IDs @ 20, 10, and 5 to prepare for persistent cluster ID assignments
all_current_20  = df.filter(pl.col("cluster_distance") == 20).select(["cluster_id", "sample_id"])
all_current_10  = df.filter(pl.col("cluster_distance") == 10).select(["cluster_id", "sample_id"])
all_current_5   = df.filter(pl.col("cluster_distance") == 5).select(["cluster_id", "sample_id"])
all_current_unclustered = df.filter(pl.col("cluster_distance") == -1).select(["cluster_id", "sample_id"])
all_current_20.write_csv('all_current_@20', sep='\t')
all_current_10.write_csv('all_current@10', sep='\t')
all_current_5.write_csv('all_current@5', sep='\t')
all_current_unclustered.write_csv('all_current@unclustered', sep='\t')




sample_clusters = renamed_samples.group_by("sample_id", maintain_order=True).agg(pl.col("cluster_id"))



exit(0)

current_samples = current_samples.with_columns([
    pl.when((pl.col("cluster_distance") == 5) & (pl.col("sample_id").is_not_null()))
      .then(current_samples.filter(pl.col("cluster_distance") == 10)
            .select(["sample_id", "cluster_id"])
            .rename({"cluster_id": "parent_cluster"}))
      .otherwise(None)
      .alias("parent_cluster")
])

current_samples = current_samples.with_columns([
    pl.when((pl.col("cluster_distance") == 10) & (pl.col("sample_id").is_not_null()))
      .then(current_samples.filter(pl.col("cluster_distance") == 20)
            .select(["sample_id", "cluster_id"])
            .rename({"cluster_id": "parent_cluster"}))
      .otherwise(None)
      .alias("parent_cluster")
])


cluster_20 = current_samples.filter(pl.col("cluster_distance") == 20)
cluster_10 = current_samples.filter(pl.col("cluster_distance") == 10)
cluster_05 = current_samples.filter(pl.col("cluster_distance") == 5)

# Write the filtered DataFrames to separate TSVs
cluster_20.write_csv(output_20_tsv, separator='\t')
cluster_10.write_csv(output_10_tsv, separator='\t')
cluster_05.write_csv(output_05_tsv, separator='\t')


exit(0)



# find common values in cluster IDs
common_values = current_ids.select(pl.col("cluster_id")).intersect(df2.select(pl.col("cluster_id")))

# Filter rows based on common values in the first column
df1_filtered = df1.filter(pl.col("column_0").is_in(common_values["column_0"]))
df2_filtered = df2.filter(pl.col("column_0").is_in(common_values["column_0"]))

# Write the filtered DataFrames to TSV files
df1_filtered.write_csv(output1_path, separator='\t', has_header=False)
df2_filtered.write_csv(output2_path, separator='\t', has_header=False)

# TODO: this is where the persistent cluster ID script needs to be

with open("temp_cluster_extraction.tsv", "w", encoding="utf-8") as temp_cluster_xtract:
    temp_cluster_xtract.write('Cluster\tSamples')
    temp_cluster_xtract.write(f"{current_cluster_id}\t{sample_ids}")
    temp_cluster_xtract.write("\n") # usher needs this or else it ignores the file
logging.info("Extracting %s subtree...", current_cluster_id)
handle_subprocess("DEBUG: temp_cluster_extraction.tsv", "cat temp_cluster_extraction.tsv")
handle_subprocess("",
    f'matUtils extract -i "{args.mat_tree}" -t "{prefix}{current_cluster_id}" -s temp_cluster_extraction.tsv -N {minimum_tree_size}')
handle_subprocess("", # TODO: restore metadata in the JSON version of the tree, -M metadata_tsv
    f'matUtils extract -i "{args.mat_tree}" -j "{prefix}{current_cluster_id}" -s temp_cluster_extraction.tsv -N {minimum_tree_size}')

# for some reason, nwk subtrees seem to end up with .nw as their extension
logging.debug("Workdir as current")
logging.debug(os.listdir('.'))
os.rename(f"{prefix}.nw", f"{prefix}.nwk")
logging.debug("Workdir, renamed nwk")
logging.debug(os.listdir('.'))

# TODO: rename persistent clusters first... but that might be a mess to do while recursing, so maybe pull this out later?
if args.microreact:
    handle_subprocess(f"Uploading {cluster_name} to MR...",
        f"python3 scripts/microreact.py {cluster_name} {prefix}_{cluster_name}.nwk {prefix}_{cluster_name}_dmtrx.tsv {samples_in_cluster_str} {args.microreacttokenfile}")
