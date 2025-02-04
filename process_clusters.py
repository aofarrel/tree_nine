import argparse
import logging
from datetime import date
import subprocess
import polars as pl # this is overkill and takes forever to import; too bad!
pl.Config.set_tbl_rows(200)


parser = argparse.ArgumentParser(description="Crunch data, extract trees, upload to MR, etc")
parser.add_argument('-ls', '--latestsamples', type=str, help='TSV: latest sample information')
#parser.add_argument('-sm', '--samplemeta', type=str, help='TSV: sample metadata pulled from terra (including myco outs), one line per sample')
#parser.add_argument('-pcm', '--persistentclustermeta', type=str, help='TSV: persistent cluster metadata from last full run of TB-D')
#parser.add_argument('-ccm', '--latestclustermeta', type=str, help='TSV: latest cluster metadata (as identified by find_clusters.py like five minutes ago)')
parser.add_argument('-pid', '--persistentids', type=str, help='TSV: persistent IDs from last full run of TB-D')
#parser.add_argument('-cid', '--latestids', type=str, help='TSV: latest cluster IDs (as identified by find_clusters.py like five minutes ago)')
#parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')

args = parser.parse_args()
all_latest_samples = pl.read_csv(args.latestsamples, separator="\t")
all_persistent_samples = pl.read_csv(args.persistentids, separator="\t")
#latest_clusters = pl.read_csv(args.latestclustermeta, separator="\t")
#persistent_clusters = pl.read_csv(args.persistentclustermeta, separator="\t")

print("all latest samples:")
print(all_latest_samples)
print("all persistent samples:")
print(all_persistent_samples)


# ensure each sample in latest-clusters has, at most, one 20 SNP, one 10 SNP, and one 05 SNP
check_clusters_valid = all_latest_samples.group_by("sample_id", maintain_order=True).agg(pl.col("cluster_distance"))
for row in check_clusters_valid.iter_rows(named=True):
    if len(row["cluster_distance"]) == 1:
        if row["cluster_distance"] == [-1]: # unclustered -- but not latestly in find_clusters.py
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

# cluster IDs @ 20, 10, and 5 to prepare for persistent cluster ID assignments
all_latest_20  = all_latest_samples.filter(pl.col("cluster_distance") == 20).select(["sample_id", "latest_cluster_id"])
all_latest_10  = all_latest_samples.filter(pl.col("cluster_distance") == 10).select(["sample_id", "latest_cluster_id"])
all_latest_5   = all_latest_samples.filter(pl.col("cluster_distance") == 5).select(["sample_id", "latest_cluster_id"])
all_latest_unclustered = all_latest_samples.filter(pl.col("cluster_distance") == -1).select(["sample_id", "latest_cluster_id"])

print("All latest 20:")
print(all_latest_20)
print("All latest 10:")
print(all_latest_10)
print("All latest 5:")
print(all_latest_5)

all_persistent_20 = all_persistent_samples.filter(pl.col("cluster_distance") == 20).select(["sample_id", "cluster_id"])
all_persistent_10  = all_persistent_samples.filter(pl.col("cluster_distance") == 10).select(["sample_id", "cluster_id"])
all_persistent_5   = all_persistent_samples.filter(pl.col("cluster_distance") == 5).select(["sample_id", "cluster_id"])
all_persistent_unclustered = all_persistent_samples.filter(pl.col("cluster_distance") == -1).select(["sample_id", "cluster_id"])

print("All persistent 20:")
print(all_persistent_20)
print("All persistent 1O:")
print(all_persistent_10)
print("All persistent 5:")
print(all_persistent_5)

# Marc's script requires that you input only sample IDs that are present in both the persistent cluster file 
# and your latest clusters, so we need to do an inner join first -- after getting our "rosetta stone" we will
# modify the original dataframe.
# "But Ash!" I hear you say, "These merges give you two dataframes that each have a column of old IDs and a
# column of new IDs! That's a rosetta stone already, we don't need Marc's script!"
# You are a fool. Yes, we could stick to that... but then we wouldn't be able to handle situations where
# clusters merge, split, or generally get messy without reinventing the wheel Marc has already made for us.
filtered_latest_20 = all_latest_20.join(all_persistent_20.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'})
filtered_latest_10 = all_latest_10.join(all_persistent_10.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'})
filtered_latest_5 = all_latest_5.join(all_persistent_5.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'})

print("filtered_latest_20:")
print(filtered_latest_20)
print("filtered_latest_10:")
print(filtered_latest_10)
print("filtered_latest_5:")
print(filtered_latest_5)

filtered_persistent_20 = all_persistent_20.join(all_latest_20.drop(['latest_cluster_id']), on="sample_id", how="inner")
filtered_persistent_10 = all_persistent_10.join(all_latest_10.drop(['latest_cluster_id']), on="sample_id", how="inner")
filtered_persistent_5 = all_persistent_5.join(all_latest_5.drop(['latest_cluster_id']), on="sample_id", how="inner")

print("filtered_persistent_20:")
print(filtered_persistent_20)
print("filtered_persistent_10:")
print(filtered_persistent_10)
print("filtered_persistent_5:")
print(filtered_persistent_5)

filtered_latest_20.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_20.tsv', separator='\t', include_header=False)
filtered_latest_10.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_10.tsv', separator='\t', include_header=False)
filtered_latest_5.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_5.tsv', separator='\t', include_header=False)

filtered_persistent_20.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_20.tsv', separator='\t', include_header=False)
filtered_persistent_10.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_10.tsv', separator='\t', include_header=False)
filtered_persistent_5.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_5.tsv', separator='\t', include_header=False)

subprocess.run("perl /scripts/marcs_incredible_script.pl filtered_persistent_20.tsv filtered_latest_20.tsv", shell=True, check=True, capture_output=True, text=True)
subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_10.tsv", shell=True, check=True)
subprocess.run("perl /scripts/marcs_incredible_script.pl filtered_persistent_10.tsv filtered_latest_10.tsv", shell=True, check=True, capture_output=True, text=True)
subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_20.tsv", shell=True, check=True)
subprocess.run("perl /scripts/marcs_incredible_script.pl filtered_persistent_5.tsv filtered_latest_5.tsv", shell=True, check=True, capture_output=True, text=True)
subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_5.tsv", shell=True, check=True)

rosetta_20 = pl.read_csv("rosetta_stone_20.tsv", separator="\t", has_header=False).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})
rosetta_10 = pl.read_csv("rosetta_stone_10.tsv", separator="\t", has_header=False).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})
rosetta_5 = pl.read_csv("rosetta_stone_20.tsv", separator="\t", has_header=False).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})

cool_samples = all_latest_samples.join(rosetta_20, on="latest_cluster_id", how="full")
cool_samples = all_latest_samples.join(rosetta_10, on="latest_cluster_id", how="full")
cool_samples = all_latest_samples.join(rosetta_5, on="latest_cluster_id", how="full")

cool_samples = cool_samples.with_columns(
    pl.when(cool_samples["sample_id"].is_in(all_persistent_20["sample_id"]))
    .then(True)
    .otherwise(False)
    .alias("in_20_cluster_last_run") # NOT AN INDICATION OF BEING BRAND NEW/NEVER CLUSTERED BEFORE
)

cool_samples = cool_samples.with_columns(
    pl.when(cool_samples["sample_id"].is_in(all_persistent_10["sample_id"]))
    .then(True)
    .otherwise(False)
    .alias("in_10_cluster_last_run") # NOT AN INDICATION OF BEING BRAND NEW/NEVER CLUSTERED BEFORE
)

cool_samples = cool_samples.with_columns(
    pl.when(cool_samples["sample_id"].is_in(all_persistent_5["sample_id"]))
    .then(True)
    .otherwise(False)
    .alias("in_5_cluster_last_run") # NOT AN INDICATION OF BEING BRAND NEW/NEVER CLUSTERED BEFORE
)


# TODO: PERSISTENT IDS ONLY WORKING FOR 10 SNP???
print("before pl.coalesce")
print(cool_samples)

cool_samples = cool_samples.with_columns(pl.coalesce('persistent_cluster_id', 'latest_cluster_id').alias("cluster_id"))
cool_samples = cool_samples.drop(['latest_cluster_id_right', 'persistent_cluster_id', 'latest_cluster_id'])

# Check for B.S.
true_for_10_not_20 = cool_samples.filter(
    pl.col("in_10_cluster_last_run") & ~pl.col("in_20_cluster_last_run")
)["sample_id"].to_list()
if true_for_10_not_20:
    raise ValueError(f"These samples were in a 10 SNP cluster last time, but not a 20 SNP cluster: {', '.join(true_for_10_not_20)}")

true_for_5_not_10 = cool_samples.filter(
    pl.col("in_5_cluster_last_run") & ~pl.col("in_10_cluster_last_run")
)["sample_id"].to_list()
if true_for_5_not_10:
    raise ValueError(f"These samples were in a 5 SNP cluster last time, but not a 10 SNP cluster: {', '.join(true_for_5_not_10)}")

true_for_5_not_20 = cool_samples.filter(
    pl.col("in_5_cluster_last_run") & ~pl.col("in_20_cluster_last_run")
)["sample_id"].to_list()
if true_for_5_not_20:
    raise ValueError(f"These samples were in a 5 SNP cluster last time, but not a 20 SNP cluster: {', '.join(true_for_5_not_20)}")

print("Cool samples")
print(cool_samples)


grouped = cool_samples.group_by("cluster_id").agg(
    pl.col("sample_id"),
    pl.col("cluster_distance").n_unique().alias("distance_nunique"),
    pl.col("cluster_distance").unique().alias("distance_values"),
    pl.col("in_20_cluster_last_run").unique(),
    pl.col("in_10_cluster_last_run").unique(),
    pl.col("in_5_cluster_last_run").unique(),
)
if (grouped["distance_nunique"] > 1).any():
    print(grouped)
    raise ValueError("Some clusters have multiple unique cluster_distance values.")

grouped = grouped.with_columns(
    grouped["distance_values"].list.get(0).alias("cluster_distance")
).drop(["distance_nunique", "distance_values"])

print("After intager-a-fy")
print(grouped)

grouped = grouped.with_columns(
    pl.when(pl.col("cluster_distance") == 20)
    .then(pl.col("in_20_cluster_last_run"))
    .otherwise(
        pl.when(pl.col("cluster_distane") == 10)
        .then(pl.col("in_10_cluster_last_run"))
        .otherwise(
            pl.when(pl.col("cluster_distance") == 5)
            .then(pl.col("in_5_cluster_last_run"))
            .otherwise(None)
        )
    )
    .alias("has_new_samples")
)

print("after looking for new samples")
print(grouped)

df = cool_samples.drop("cluster_distance").join(grouped, on="cluster_id")
print("After join")
print(df)
df = df.with_columns(
    pl.lit(None).cast(pl.Utf8).alias("parent"),
    pl.lit([]).cast(pl.List(pl.Utf8)).alias("children")
)
df = df.sort("cluster_distance")
print("sorted")
print(df)

    
# Mapping from sample_id lists to cluster_id by distance
sample_map = {dist: {} for dist in [5, 10, 20]}
for row in df.iter_rows(named=True):
    sample_map[row["cluster_distance"]][frozenset(row["sample_id"])] = row["cluster_id"]
print("sample_map")
print(sample_map)
updates = []
for row in df.iter_rows(named=True):
    cluster_id, samples, distance = row["cluster_id"], frozenset(row["sample_id"]), row["cluster_distance"]
    
    if distance == 5:
        parent_id = sample_map[10].get(samples)
        if parent_id:
            updates.append((cluster_id, "parent", parent_id))
            updates.append((parent_id, "children", cluster_id))
    elif distance == 10:
        parent_id = sample_map[20].get(samples)
        if parent_id:
            updates.append((cluster_id, "parent", parent_id))
            updates.append((parent_id, "children", cluster_id))

# Apply updates
for cluster_id, col, value in updates:
    if col == "parent":
        df = df.with_columns(
            pl.when(df["cluster_id"] == cluster_id).then(value).otherwise(df["parent"]).alias("parent")
        )
    else:
        df = df.with_columns(
            pl.when(df["cluster_id"] == cluster_id).then(df["children"] + [value]).otherwise(df["children"]).alias("children")
        )
print(df)


#all_latest_20.write_csv('all_latest_@20', separator='\t')
#all_latest_10.write_csv('all_latest@10', separator='\t')
#all_latest_5.write_csv('all_latest@5', separator='\t')
#all_latest_unclustered.write_csv('all_latest@unclustered', separator='\t')




#sample_clusters = renamed_samples.group_by("sample_id", maintain_order=True).agg(pl.col("cluster_id"))




'''
latest_samples = latest_samples.with_columns([
    pl.when((pl.col("cluster_distance") == 5) & (pl.col("sample_id").is_not_null()))
      .then(latest_samples.filter(pl.col("cluster_distance") == 10)
            .select(["sample_id", "cluster_id"])
            .rename({"cluster_id": "parent_cluster"}))
      .otherwise(None)
      .alias("parent_cluster")
])

latest_samples = latest_samples.with_columns([
    pl.when((pl.col("cluster_distance") == 10) & (pl.col("sample_id").is_not_null()))
      .then(latest_samples.filter(pl.col("cluster_distance") == 20)
            .select(["sample_id", "cluster_id"])
            .rename({"cluster_id": "parent_cluster"}))
      .otherwise(None)
      .alias("parent_cluster")
])


cluster_20 = latest_samples.filter(pl.col("cluster_distance") == 20)
cluster_10 = latest_samples.filter(pl.col("cluster_distance") == 10)
cluster_05 = latest_samples.filter(pl.col("cluster_distance") == 5)

# Write the filtered DataFrames to separate TSVs
cluster_20.write_csv(output_20_tsv, separator='\t')
cluster_10.write_csv(output_10_tsv, separator='\t')
cluster_05.write_csv(output_05_tsv, separator='\t')


exit(0)



# find common values in cluster IDs
common_values = latest_ids.select(pl.col("cluster_id")).intersect(df2.select(pl.col("cluster_id")))

# Filter rows based on common values in the first column
df1_filtered = df1.filter(pl.col("column_0").is_in(common_values["column_0"]))
df2_filtered = df2.filter(pl.col("column_0").is_in(common_values["column_0"]))

# Write the filtered DataFrames to TSV files
df1_filtered.write_csv(output1_path, separator='\t', has_header=False)
df2_filtered.write_csv(output2_path, separator='\t', has_header=False)

# TODO: this is where the persistent cluster ID script needs to be

with open("temp_cluster_extraction.tsv", "w", encoding="utf-8") as temp_cluster_xtract:
    temp_cluster_xtract.write('Cluster\tSamples')
    temp_cluster_xtract.write(f"{latest_cluster_id}\t{sample_ids}")
    temp_cluster_xtract.write("\n") # usher needs this or else it ignores the file
logging.info("Extracting %s subtree...", latest_cluster_id)
handle_subprocess("DEBUG: temp_cluster_extraction.tsv", "cat temp_cluster_extraction.tsv")
handle_subprocess("",
    f'matUtils extract -i "{args.mat_tree}" -t "{prefix}{latest_cluster_id}" -s temp_cluster_extraction.tsv -N {minimum_tree_size}')
handle_subprocess("", # TODO: restore metadata in the JSON version of the tree, -M metadata_tsv
    f'matUtils extract -i "{args.mat_tree}" -j "{prefix}{latest_cluster_id}" -s temp_cluster_extraction.tsv -N {minimum_tree_size}')

# for some reason, nwk subtrees seem to end up with .nw as their extension
logging.debug("Workdir as latest")
logging.debug(os.listdir('.'))
os.rename(f"{prefix}.nw", f"{prefix}.nwk")
logging.debug("Workdir, renamed nwk")
logging.debug(os.listdir('.'))

# TODO: rename persistent clusters first... but that might be a mess to do while recursing, so maybe pull this out later?
if args.microreact:
    handle_subprocess(f"Uploading {cluster_name} to MR...",
        f"python3 scripts/microreact.py {cluster_name} {prefix}_{cluster_name}.nwk {prefix}_{cluster_name}_dmtrx.tsv {samples_in_cluster_str} {args.microreacttokenfile}")
'''