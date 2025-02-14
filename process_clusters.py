VERSION = "0.0.1"
verbose = True
print(f"PROCESS CLUSTERS - VERSION {VERSION}")

# pylint: disable=W0311,W1514,C0103,C0321,C0301,C0413

# Note to future maintainers: This script creates a dataframe that is extremely redundant. Worst case, every
# sample will be present three times -- once per cluster. Additionally, we are storing cluster information in
# this dataframe too, so every sample in a cluster has duplicate information about said cluster. This is
# laughably inefficient, but it might not be worth refactoring. We are using "polars" here, which is like pandas
# but significantly more efficient. Based on my experience working with Literally Every Mycobacterium Sample On
# NCBI SRA's And Its Metadata, I estimate we can get away with this ineffecient method I have written here unless
# we end up with well over 50,000 samples.

import io
import os
import csv
import json
import argparse
from datetime import datetime
import subprocess
import requests
import polars as pl # this is overkill and takes forever to import; too bad!
pl.Config.set_tbl_rows(160)
pl.Config.set_tbl_cols(-1)
pl.Config.set_tbl_width_chars(200)
pl.Config.set_fmt_str_lengths(50)
pl.Config.set_fmt_table_cell_list_len(5)
today = datetime.utcnow().date().isoformat() # I don't care if this runs past midnight, give everything the same day!
print(today)

parser = argparse.ArgumentParser(description="Crunch data, extract trees, upload to MR, etc")
parser.add_argument('-to', '--token', type=str, required=True, help="TXT: MR token")
parser.add_argument('-ls', '--latestsamples', type=str, help='TSV: latest sample information')
#parser.add_argument('-sm', '--samplemeta', type=str, help='TSV: sample metadata pulled from terra (including myco outs), one line per sample')
parser.add_argument('-pcm', '--persistentclustermeta', type=str, help='TSV: persistent cluster metadata from last full run of TB-D')
#parser.add_argument('-ccm', '--latestclustermeta', type=str, help='TSV: latest cluster metadata (as identified by find_clusters.py like five minutes ago)')
parser.add_argument('-pid', '--persistentids', type=str, help='TSV: persistent IDs from last full run of TB-D')
#parser.add_argument('-cid', '--latestids', type=str, help='TSV: latest cluster IDs (as identified by find_clusters.py like five minutes ago)')
#parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')
parser.add_argument('-mat', '--mat_tree', type=str, help='PB: tree')
parser.add_argument('-cs', '--contextsamples', type=int, default=0, help="int: Number of context samples for cluster subtrees")

args = parser.parse_args()
all_latest_samples = pl.read_csv(args.latestsamples, separator="\t")
all_persistent_samples = pl.read_csv(args.persistentids, separator="\t")
persistent_clusters_meta = pl.read_csv(args.persistentclustermeta, separator="\t", null_values="NULL", try_parse_dates=True)
#latest_clusters = pl.read_csv(args.latestclustermeta, separator="\t")
with open(args.token, 'r') as file:
    token = file.readline()

if verbose:
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

if verbose:
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

if verbose:
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

if verbose:
    print("filtered_latest_20:")
    print(filtered_latest_20)
    print("filtered_latest_10:")
    print(filtered_latest_10)
    print("filtered_latest_5:")
    print(filtered_latest_5)

filtered_persistent_20 = all_persistent_20.join(all_latest_20.drop(['latest_cluster_id']), on="sample_id", how="inner")
filtered_persistent_10 = all_persistent_10.join(all_latest_10.drop(['latest_cluster_id']), on="sample_id", how="inner")
filtered_persistent_5 = all_persistent_5.join(all_latest_5.drop(['latest_cluster_id']), on="sample_id", how="inner")

if verbose:
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
subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_20.tsv", shell=True, check=True)
subprocess.run("perl /scripts/marcs_incredible_script.pl filtered_persistent_10.tsv filtered_latest_10.tsv", shell=True, check=True, capture_output=True, text=True)
subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_10.tsv", shell=True, check=True)
subprocess.run("perl /scripts/marcs_incredible_script.pl filtered_persistent_5.tsv filtered_latest_5.tsv", shell=True, check=True, capture_output=True, text=True)
subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_5.tsv", shell=True, check=True)

rosetta_20 = pl.read_csv("rosetta_stone_20.tsv", separator="\t", has_header=False).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})
rosetta_10 = pl.read_csv("rosetta_stone_10.tsv", separator="\t", has_header=False).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})
rosetta_5 = pl.read_csv("rosetta_stone_5.tsv", separator="\t", has_header=False).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})

latest_samples_translated = (all_latest_samples.join(rosetta_20, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_20_cluster_id'}).drop("latest_cluster_id_right")
latest_samples_translated = (latest_samples_translated.join(rosetta_10, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_10_cluster_id'}).drop("latest_cluster_id_right")
latest_samples_translated = (latest_samples_translated.join(rosetta_5, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_5_cluster_id'}).drop("latest_cluster_id_right")
all_latest_samples = None

latest_samples_translated = latest_samples_translated.with_columns(
    pl.when(pl.col("cluster_distance") == 20)
    .then(
        pl.when(latest_samples_translated["sample_id"].is_in(all_persistent_20["sample_id"]))
        .then(True)
        .otherwise(False)
    )
    .otherwise(None)
    .alias("in_20_cluster_last_run") # NOT AN INDICATION OF BEING BRAND NEW/NEVER CLUSTERED BEFORE
)

latest_samples_translated = latest_samples_translated.with_columns(
    pl.when(pl.col("cluster_distance") == 10)
    .then(
        pl.when(latest_samples_translated["sample_id"].is_in(all_persistent_10["sample_id"]))
        .then(True)
        .otherwise(False)
    )
    .otherwise(None)
    .alias("in_10_cluster_last_run") # NOT AN INDICATION OF BEING BRAND NEW/NEVER CLUSTERED BEFORE
)

latest_samples_translated = latest_samples_translated.with_columns(
    pl.when(pl.col("cluster_distance") == 5)
    .then(
        pl.when(latest_samples_translated["sample_id"].is_in(all_persistent_5["sample_id"]))
        .then(True)
        .otherwise(False)
    )
    .otherwise(None)
    .alias("in_5_cluster_last_run") # NOT AN INDICATION OF BEING BRAND NEW/NEVER CLUSTERED BEFORE
)

if verbose:
    print("Before pl.coalesce")
    print(latest_samples_translated)

latest_samples_translated = latest_samples_translated.with_columns(
    pl.coalesce('persistent_20_cluster_id', 'persistent_10_cluster_id', 'persistent_5_cluster_id', 'latest_cluster_id')
    .alias("cluster_id")
)
latest_samples_translated = latest_samples_translated.drop(['persistent_20_cluster_id', 'persistent_10_cluster_id', 'persistent_5_cluster_id'])
latest_samples_translated = latest_samples_translated.rename({'latest_cluster_id': 'workdir_cluster_id'})

# Check for B.S.
true_for_10_not_20 = latest_samples_translated.filter(
    pl.col("in_10_cluster_last_run") & ~pl.col("in_20_cluster_last_run")
)["sample_id"].to_list()
if true_for_10_not_20:
    raise ValueError(f"These samples were in a 10 SNP cluster last time, but not a 20 SNP cluster: {', '.join(true_for_10_not_20)}")

true_for_5_not_10 = latest_samples_translated.filter(
    pl.col("in_5_cluster_last_run") & ~pl.col("in_10_cluster_last_run")
)["sample_id"].to_list()
if true_for_5_not_10:
    raise ValueError(f"These samples were in a 5 SNP cluster last time, but not a 10 SNP cluster: {', '.join(true_for_5_not_10)}")

true_for_5_not_20 = latest_samples_translated.filter(
    pl.col("in_5_cluster_last_run") & ~pl.col("in_20_cluster_last_run")
)["sample_id"].to_list()
if true_for_5_not_20:
    raise ValueError(f"These samples were in a 5 SNP cluster last time, but not a 20 SNP cluster: {', '.join(true_for_5_not_20)}")

if verbose:
    print("Cool samples")
    print(latest_samples_translated)

# group by persistent cluster ID
# TODO: How does this affect unclustered samples?
grouped = latest_samples_translated.group_by("cluster_id").agg(
    pl.col("sample_id"),
    pl.col("cluster_distance").n_unique().alias("distance_nunique"),
    pl.col("cluster_distance").unique().alias("distance_values"),
    pl.col("in_20_cluster_last_run").unique(),
    pl.col("in_10_cluster_last_run").unique(),
    pl.col("in_5_cluster_last_run").unique(),
    pl.col("workdir_cluster_id").unique(),
)
if (grouped["distance_nunique"] > 1).any():
    print(grouped)
    raise ValueError("Some clusters have multiple unique cluster_distance values.")

grouped = grouped.with_columns(
    grouped["distance_values"].list.get(0).alias("cluster_distance")
).drop(["distance_nunique", "distance_values"])

if verbose:
    print("After grouping and then intager-a-fy")
    print(grouped)

grouped = grouped.with_columns(
    pl.when(pl.col("cluster_distance") == 20)
    .then(pl.col("in_20_cluster_last_run"))
    .otherwise(
        pl.when(pl.col("cluster_distance") == 10)
        .then(pl.col("in_10_cluster_last_run"))
        .otherwise(
            pl.when(pl.col("cluster_distance") == 5)
            .then(pl.col("in_5_cluster_last_run"))
            .otherwise(None)
        )
    )
    .alias("samples_previously_in_cluster")
)

if verbose:
    # drop AFTER this print so we can see if anything is whacky here
    print("After looking for new samples")
    print(grouped)

grouped = grouped.drop(['in_20_cluster_last_run', 'in_10_cluster_last_run', 'in_5_cluster_last_run']) # will be readded upon join
grouped = grouped.drop("sample_id") # prevent creation of sample_id_right, also this is redundant when we agg() again later
hella_redundant = (latest_samples_translated.drop("cluster_distance")).join(grouped, on="cluster_id")
grouped = None
hella_redundant = hella_redundant.with_columns(
    pl.lit(None).cast(pl.Utf8).alias("cluster_parent"),
    pl.lit([]).cast(pl.List(pl.Utf8)).alias("cluster_children")
)
hella_redundant = hella_redundant.sort("cluster_distance")
latest_samples_translated = None

# we will be using hella_redundant again later
if verbose:
    print("Joined grouped with latest_samples_translated to form hella_redundant, and sorted by cluster distance")
    print(hella_redundant)

# this is a super goofy way to link parents and children, but it seems to work
# basically, we're building a dataframe that has per-sample information, and a bunch of
# duplicated per-cluster information. this is EXTREMELY redundant and in a normal world
# we would not do this!
sample_map = {dist: {} for dist in [5, 10, 20]}
for row in hella_redundant.iter_rows(named=True):
    sample_map[row["cluster_distance"]][row["sample_id"]] = row["cluster_id"]
if verbose:
    print("Sample map")
    print(sample_map)
updates = []
for row in hella_redundant.iter_rows(named=True):
    cluster_id, one_sample, distance = row["cluster_id"], row["sample_id"], row["cluster_distance"]
    print(f"[{distance}] sample {one_sample} in cluster_id {cluster_id}")
    if distance == 5:
        parent_id = sample_map[10].get(one_sample)
        print(f"parent_id {parent_id}")
        if parent_id:
            updates.append((cluster_id, "cluster_parent", parent_id))
    elif distance == 10:
        parent_id = sample_map[20].get(one_sample)
        print(f"parent_id {parent_id}")
        if parent_id:
            updates.append((cluster_id, "cluster_parent", parent_id))
        child_id = sample_map[5].get(one_sample)
        print(f"child_id {child_id}") # ONLY GRABS THE CHILD ID ASSOCIATED WITH THIS SAMPLE, NOT ALL CHILDREN OF CLUSTER
        if child_id:
            updates.append((cluster_id, "cluster_one_child", child_id))
    elif distance == 20:
        child_id = sample_map[10].get(one_sample)
        print(f"child_id {child_id}") # ONLY GRABS THE CHILD ID ASSOCIATED WITH THIS SAMPLE, NOT ALL CHILDREN OF CLUSTER
        if child_id:
            updates.append((cluster_id, "cluster_one_child", child_id))
    else:
        raise ValueError
if verbose:
    print("updates")
    print(updates)

for cluster_id, col, value in updates:
    print(f"For cluster {cluster_id}, col {col}, val {value} in updates")
    if col == "cluster_parent":
        hella_redundant = hella_redundant.with_columns(
            pl.when(pl.col("cluster_id") == cluster_id)
            .then(pl.lit(value))
            .otherwise(pl.col("cluster_parent"))
            .alias("cluster_parent")
        )
    else:
        hella_redundant = hella_redundant.with_columns(
            pl.when(pl.col("cluster_id") == cluster_id)
            .then((pl.col("cluster_children").list.concat(pl.lit(value))).list.unique())
            .otherwise(pl.col("cluster_children"))
            .alias("cluster_children")
        )
cluster_id = None
if verbose:
    print("after linking parents and children")
    print(hella_redundant)

hella_redundant = hella_redundant.with_columns([
    # When samples_previously_in_cluster is:
    # [True, False]/[False, True] --> some samples were in cluster previously     --> old cluster, needs updating
    # [False]                     --> no samples were in this cluster previously  --> new cluster, needs updating
    # [True]                      --> all samples were in this cluster previously --> old cluster, unchanged
    #
    # Strictly speaking you don't need to update something that didn't exist anymore but I gave this a good long think
    # and have decided that "needs updating" (ie pinging microreact) is functionally what matters here.
    pl.when(pl.col("samples_previously_in_cluster") == [False])
    .then(True)
    .otherwise(False)
    .alias("cluster_brand_new"),

    pl.when(pl.col("samples_previously_in_cluster") == [True])
    .then(False)
    .otherwise(True)
    .alias("cluster_needs_updating"),

    # this particcular row is PER SAMPLE, not PER ClUSTER
    ~pl.coalesce(["in_20_cluster_last_run", "in_10_cluster_last_run", "in_5_cluster_last_run"]).alias("sample_newly_clustered")
]).drop("samples_previously_in_cluster")

if verbose:
    print("after processing what clusters are brand new and what ones are unchanged")
    print(hella_redundant)

# drop AFTER the verbose print
hella_redundant = hella_redundant.drop(["in_20_cluster_last_run", "in_10_cluster_last_run", "in_5_cluster_last_run"])

# We're grouping again! Is there a way to do this in the previous group? Maybe, but I'm trying
# to get this up and running ASAP so whatever works, works!
second_group = hella_redundant.group_by("cluster_id").agg(
    pl.col("sample_id"),
    pl.col("cluster_distance").unique(),
    pl.col("cluster_brand_new").unique(),
    pl.col("cluster_needs_updating").unique(),
    pl.col("cluster_parent").unique(),
    pl.col("cluster_children").flatten().unique(),
    pl.col("workdir_cluster_id").unique().first(),
)
if verbose:
    print("after grouping hella_redundant by cluster_id")
    print(second_group)

# check cluster distances
# TODO: this and other asserts will probably need to change if we change how we handle unclustered samples
assert (second_group["cluster_distance"].list.len() == 1).all(), "cluster_distance lists have length ≠ 1"
second_group = second_group.with_columns(pl.col("cluster_distance").list.get(0).alias("cluster_distance_int"))
second_group = second_group.drop("cluster_distance").rename({"cluster_distance_int": "cluster_distance"})
print(second_group)

# check parenthood
# in polars, [null] is considered to have a length of 1. I'm checking by distance rather than all clusters at once in case that changes.
# be aware: https://github.com/pola-rs/polars/issues/18522
assert ((second_group.filter(second_group["cluster_distance"] == 5))["cluster_parent"].list.len() == 1).all(), "5-cluster with multiple cluster_parents"
assert ((second_group.filter(second_group["cluster_distance"] == 10))["cluster_parent"].list.len() == 1).all(), "10-cluster with multiple cluster_parents"
assert ((second_group.filter(second_group["cluster_distance"] == 20))["cluster_parent"].list.len() == 1).all(), "20-cluster with cluster_parent *OR* preliminary null error"
second_group = second_group.with_columns(pl.col("cluster_parent").list.get(0).alias("cluster_parent_str"))
second_group = second_group.drop("cluster_parent").rename({"cluster_parent_str": "cluster_parent"})
assert ((second_group.filter(second_group["cluster_distance"] == 20))["cluster_parent"].is_null()).all(), "20-cluster with cluster_parent *OR* secondary null error"

# check otherstuff
assert ((second_group.filter(second_group["cluster_distance"] == 5))["cluster_children"].list.get(0).is_null()).all(), "5-cluster with cluster_children"
assert (second_group["cluster_needs_updating"].list.len() == 1).all(), "Cluster not sure if it needs updating"
second_group = second_group.with_columns(pl.col("cluster_needs_updating").list.get(0).alias("cluster_needs_updating_bool"))
second_group = second_group.drop("cluster_needs_updating").rename({"cluster_needs_updating_bool": "cluster_needs_updating"})
assert (second_group["cluster_brand_new"].list.len() == 1).all(), "Cluster not sure if it's new"
second_group = second_group.with_columns(pl.col("cluster_brand_new").list.get(0).alias("cluster_brand_new_bool"))
second_group = second_group.drop("cluster_brand_new").rename({"cluster_brand_new_bool": "cluster_brand_new"})

# convert [null] to null
second_group = second_group.with_columns(
    pl.when(pl.col("cluster_children").list.get(0).is_null())
    .then(None)
    .otherwise(pl.col("cluster_children"))
    .alias("cluster_children")
)

# join with persistent cluster metadata tsv
# TODO: eventually latest cluster metadata file should be joined here too
all_cluster_information = second_group.join(persistent_clusters_meta, how="full", on="cluster_id", coalesce=True)

# TODO: this would be cool but also you could just sort the data table by the date-collected or date-added column probably?
# now let's get information as to which samples are new or old so we can highlight them\
#for row in hella_redundant.iter_rows(named=True):
#    this_cluster_id = row["cluster_id"]
#    this_sample_id = row["sample_id"]
#    distance = row["cluster_distance"]
#    in_20 = row["in_20_cluster_last_run"]
#    in_10 = row["in_10_cluster_last_run"]
#    in_5 = row["in_5_cluster_last_run"]
#    brand_new_sample = True if this_sample_id in all_persistent_samples["sample_id"] else False

if "last_update" not in all_cluster_information.columns:
    print("No last_update column, adding an empty one...")
    all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("last_update"))
if "first_found" not in all_cluster_information.columns:
    print("No first_found column, adding an empty one...")
    all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("first_found"))

if verbose:
    print("Old persistent cluster metadata")
    print(persistent_clusters_meta)
    print("All cluster information we have so far after joining with that")
    print(all_cluster_information)

# okay, everything looks good so far. let's get some URLs!!
print("Assigning self-URLs...")
for row in all_cluster_information.iter_rows(named=True):
    this_cluster_id = row["cluster_id"]
    distance = row["cluster_distance"]
    has_children = False if row["cluster_children"] is None else True
    brand_new = row["cluster_brand_new"]
    needs_updating = row["cluster_needs_updating"]
    URL = row["microreact_url"]

    # this is the only type of cluster that doesn't need a URL
    # TODO: currently, as written, this means we don't really keep track of 20-no-kids that gain/lose samples! is this good?
    if distance == 20 and not has_children:
        print(f"{this_cluster_id} is 20-cluster with no children (but maybe new samples), skipping")
        continue

    if brand_new:
        assert URL is None, f"{this_cluster_id} is brand new but already has a MR URL?"
        first_found = today
        last_update = today

        with open("./BLANK_template.json", "r") as temp_proj_json:
            mr_document = json.load(temp_proj_json)
        update_resp = requests.post("https://microreact.org/api/projects/create",
            headers={"Access-Token": token, "Content-Type": "application/json; charset=UTF-8"},
            params={"access": "private"},
            timeout=100,
            json=mr_document)
        if update_resp.status_code == 200:
            URL = update_resp.json()['id']
            if verbose: print(f"Created temporary MR project with id {URL}")
        else:
            print("Failed to create new project:", update_resp.status_code, update_resp.text)
        print(f"{this_cluster_id} is brand new and has been assigned MR URL {URL} and a first-found date")
        all_cluster_information = all_cluster_information.with_columns(
            pl.when(all_cluster_information["cluster_id"] == this_cluster_id)
            .then(pl.lit(first_found))
            .otherwise(all_cluster_information["first_found"])
            .alias("first_found"))

    elif needs_updating:
        # later on, assert URL and first_found is not None... but for the first time don't do that!
        if URL is None:
            print(f"WARNING: {this_cluster_id} isn't brand new, but is flagged as needing an update and has no URL. Will make a new URL.")
            with open("./BLANK_template.json", "r") as temp_proj_json:
                mr_document = json.load(temp_proj_json)
            update_resp = requests.post("https://microreact.org/api/projects/create",
                headers={"Access-Token": token, "Content-Type": "application/json; charset=UTF-8"},
                params={"access": "private"},
                timeout=100,
                json=mr_document)
            if update_resp.status_code == 200:
                URL = update_resp.json()['id']
                if verbose: print(f"Created temporary MR project with id {URL}")
            else:
                print("Failed to create new project:", update_resp.status_code, update_resp.text)
            print(f"{this_cluster_id} is brand new and has been assigned MR URL {URL} and a first-found date")
        else:
            # TODO: check cluster URL is valid here... unless it's gonna be a pain in the neck then maybe dont
            print(f"{this_cluster_id}'s URL is valid ({URL})")
        last_update = today

    else:
        # unchanged
        continue

    all_cluster_information = all_cluster_information.with_columns([
        pl.when(all_cluster_information["cluster_id"] == this_cluster_id)
        .then(pl.lit(URL))
        .otherwise(all_cluster_information["microreact_url"])
        .alias("microreact_url"),
        pl.when(all_cluster_information["cluster_id"] == this_cluster_id)
        .then(pl.lit(last_update))
        .otherwise(all_cluster_information["last_update"])
        .alias("last_update")
    ])

all_cluster_information = all_cluster_information.with_columns(
    pl.lit(None).alias("parent_URL"),
    pl.lit(None).alias("children_URLs") # TODO: how are gonna keep these ordered with the children column... does that even matter?
)

# now that everything has a URL, or doesn't need one, iterate a second time to get URLs of parents and children
print("Searching for MR URLs of parents and children...")
for row in all_cluster_information.iter_rows(named=True):
    this_cluster_id = row["cluster_id"]
    distance = row["cluster_distance"]
    has_children = False if row["cluster_children"] is None else True
    has_parent = False if row["cluster_parent"] is None else True
    needs_updating = row["cluster_needs_updating"]
    brand_new = row["cluster_brand_new"]
    URL = row["microreact_url"]

    # Because there is never a situation where a new child cluster pops up in a parent cluster that doesn't need to be updated,
    # and because MR URLs don't need to be updated, clusters that don't need updating don't need to know parent/child URLs.
    if not needs_updating:
        print(f"{this_cluster_id} is a {distance}-cluster with no new samples (ergo no new children), skipping")
        continue

    if has_parent:
        # get value of column "microreact_url" when column "cluster_id" matches the string in row["cluster_parent"]
        parent_URL = all_cluster_information.filter(
            all_cluster_information["cluster_id"] == row["cluster_parent"]
        )["microreact_url"].item()
        print(f"{this_cluster_id} has parent URL {parent_URL}")
        all_cluster_information = all_cluster_information.with_columns(
            pl.when(all_cluster_information["cluster_id"] == this_cluster_id)
            .then(pl.lit(parent_URL))
            .otherwise(all_cluster_information["parent_URL"])
            .alias("parent_URL")
        )
    if has_children:
        # get value of column "microreact_url" when column "cluster_id" matches a string in the list in row["cluster_children"]
        children_URLs = all_cluster_information.filter(
            all_cluster_information["cluster_id"].is_in(row["cluster_children"])
        )["microreact_url"].to_list()
        print(f"{this_cluster_id} has children with URLs {children_URLs}")
        all_cluster_information = all_cluster_information.with_columns(
            pl.when(all_cluster_information["cluster_id"] == this_cluster_id)
            .then(pl.lit(children_URLs))
            .otherwise(all_cluster_information["children_URLs"])
            .alias("children_URLs")
        )
parent_URL, children_URLs, URL = None, None, None
print("Final data table:")
print(all_cluster_information)
all_cluster_information.write_ndjson('all_cluster_information_before_upload.json')

# now that everything can be crosslinked, iterate one more time to actually upload
# yes three iterations is weird, cringe even. I'm doing this to make debugging easier.
print("Splitting and uploading...")
for row in all_cluster_information.iter_rows(named=True):
    this_cluster_id = row["cluster_id"]
    distance = row["cluster_distance"]
    sample_id_list = row["sample_id"]
    has_children = False if row["cluster_children"] is None else True
    has_parent = False if row["cluster_parent"] is None else True
    parent_URL = row["parent_URL"] # children URLs handled differently
    cluster_parent = row["cluster_parent"] # None if !has_parent
    n_children = len(row["cluster_children"]) if has_children else -1
    needs_updating = row["cluster_needs_updating"]
    needs_updating = row["cluster_needs_updating"]
    URL = row["microreact_url"]
    workdir_cluster_id = row["workdir_cluster_id"]

    # Because there is never a situation where a new child cluster pops up in a parent cluster that doesn't need to be updated,
    # and because MR URLs don't need to be updated, clusters that don't need updating don't need to know parent/child URLs.
    if not needs_updating:
        print(f"{this_cluster_id} is a {distance}-cluster with no new samples, skipping")
        continue
    
    # check for cluster-level distance matrix and nwk
    if workdir_cluster_id is None:
        with open(f"a{this_cluster_id}_dmtrx.tsv", "r") as distance_matrix:
            this_a_matrix = distance_matrix.readlines()
        with open(f"a{this_cluster_id}.nwk", "r") as nwk_file:
            this_a_nwk = nwk_file.readline() # only need first line
        
        # backmasked
        if os.path.exists(f"b{this_cluster_id}_dmtrx.tsv"):
            with open(f"b{this_cluster_id}_dmtrx.tsv", "r") as distance_matrix:
                this_b_matrix = distance_matrix.readlines()
            with open(f"b{this_cluster_id}.nwk", "r") as nwk_file:
                this_b_nwk = nwk_file.readline() # only need first line
        else:
            has_backmask = False
            this_b_matrix = "Backmasked nwk not available\tSorry\nBackmasked nwk not available\t0\t0\nSorry\t0\t0"
            this_b_nwk = '("Backmasked nwk not available","Sorry"):0;'
    
    else:
        with open(f"a{workdir_cluster_id}_dmtrx.tsv", "r") as distance_matrix:
            this_a_matrix = distance_matrix.readlines()
        with open(f"a{workdir_cluster_id}.nwk", "r") as nwk_file:
            this_a_nwk = nwk_file.readline() # only need first line

        # rename to persistent IDs so WDL-globbed outs actually makes sense
        if os.path.exists(f"a{this_cluster_id}_dmtrx.tsv"): # this should never happen and indicates an issue with persistent cluster ID
            print(f"WARNING: Successfully loaded a{workdir_cluster_id}_dmtrx.tsv but could not rename it with persistent ID a{this_cluster_id} as a dmatrix with that name already exists")
        else:
            os.rename(f"a{workdir_cluster_id}_dmtrx.tsv", f"a{this_cluster_id}_dmtrx.tsv")
            print(f"Renamed a{workdir_cluster_id}_dmtrx.tsv to its persistent ID form a{this_cluster_id}_dmtrx.tsv")
        if os.path.exists(f"a{this_cluster_id}.nwk"): # this should never happen and indicates an issue with persistent cluster ID
            print(f"WARNING: Successfully loaded a{workdir_cluster_id}.nwk but could not rename it with persistent ID a{this_cluster_id} as a nwk with that name already exists")
        else:
            os.rename(f"a{workdir_cluster_id}.nwk", f"a{this_cluster_id}.nwk")
            print(f"Renamed a{workdir_cluster_id}.nwk to its persistent ID form a{this_cluster_id}.nwk")
        
        # backmasked
        if os.path.exists(f"b{workdir_cluster_id}_dmtrx.tsv"):
            has_backmask = True
            with open(f"b{workdir_cluster_id}_dmtrx.tsv", "r") as distance_matrix:
                this_b_matrix = distance_matrix.readlines()
            with open(f"b{workdir_cluster_id}.nwk", "r") as nwk_file:
                this_b_nwk = nwk_file.readline() # only need first line

            if os.path.exists(f"b{this_cluster_id}_dmtrx.tsv"): # this should never happen and indicates an issue with persistent cluster ID
                print(f"WARNING: Successfully loaded b{workdir_cluster_id}_dmtrx.tsv but could not rename it with persistent ID b{this_cluster_id} as a dmatrix with that name already exists")
            else:
                os.rename(f"b{workdir_cluster_id}_dmtrx.tsv", f"b{this_cluster_id}_dmtrx.tsv")
                print(f"Renamed b{workdir_cluster_id}_dmtrx.tsv to its persistent ID form b{this_cluster_id}_dmtrx.tsv")
            if os.path.exists(f"b{this_cluster_id}.nwk"): # this should never happen and indicates an issue with persistent cluster ID
                print(f"WARNING: Successfully loaded b{workdir_cluster_id}.nwk but could not rename it with persistent ID b{this_cluster_id} as a nwk with that name already exists")
            else:
                os.rename(f"b{workdir_cluster_id}.nwk", f"b{this_cluster_id}.nwk")
                print(f"Renamed b{workdir_cluster_id}.nwk to its persistent ID form b{this_cluster_id}.nwk")
        else:
            has_backmask = False
            this_b_matrix = "Backmasked nwk not available\tSorry\nBackmasked nwk not available\t0\t0\nSorry\t0\t0"
            this_b_nwk = '("Backmasked nwk not available","Sorry"):0;'

    with open("./REALER_template.json", "r") as real_template_json:
        mr_document = json.load(real_template_json)

    # project title
    mr_document["meta"]["name"] = f"{this_cluster_id} updated {today}"
    mr_document["meta"]["description"] = f"{this_cluster_id} as automatically generated by version {VERSION} of process_clusters.py"

    # note
    markdown_note = f"### {this_cluster_id} ({distance}-SNP, {len(sample_id_list)} samples)\nUpdated {today}\n\n"
    if len(sample_id_list) == 2:
        markdown_note += "*WARNING: This is an extremely small cluster and its tree may not render in Microreact!\n"
    #markdown_note += "The default view shows the tree and distance matrix before intra-cluster backmasking. "
    if has_backmask:
        markdown_note += "Please click the 'backmask' tab in tree/matrix to view the backmasked versions of this cluster.\n\n"
    else:
        markdown_note += "Masking within clusters ('backmasking') is being revamped; SNP distance may be slightly inflated in the meantime.\n\n"
    if has_parent:
        markdown_note += f"Parent cluster: [{cluster_parent}](https://microreact.org/project/{parent_URL})\n\n"
    if has_children:
        markdown_note += f"Child clusters ({n_children} total):\n" + "".join(f"* [click here](https://microreact.org/project/{child_url})\n" for child_url in row["children_URLs"]) + "\n"
    mr_document["notes"]["note-1"]["source"] = markdown_note

    # trees
    mr_document["files"]["chya"]["name"] = f"a{this_cluster_id}.nwk"
    mr_document["files"]["chya"]["blob"] = this_a_nwk
    mr_document["files"]["bmtr"]["name"] = f"b{this_cluster_id}.nwk"
    mr_document["files"]["bmtr"]["blob"] = this_a_nwk # TODO: replace with this_b_nwk!!

    # tree labels
    # TODO: MR needs all panels filled out or else the workspace won't load. We're skipping metadata for now, so we're just doing not
    # very useful stuff for now. also using io like this is a little whack but it'll do
    crappy_temporary_metadata_table = [
        {
            "id": sample_id,
            "state": "probably California",
            "country": "likely USA",
        }
        for sample_id in sample_id_list
    ]
    output = io.StringIO(newline='') # get rid of carriage return
    writer = csv.DictWriter(output, fieldnames=crappy_temporary_metadata_table[0].keys(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(crappy_temporary_metadata_table)
    labels = output.getvalue()
    if verbose: print(labels)
    mr_document["files"]["ji0o"]["name"] = f"{this_cluster_id}_metadata.csv"
    mr_document["files"]["ji0o"]["blob"] = labels

    # distance matrix already readlines()'d into this_a_matrix, just make it a string
    this_a_matrix = "\n".join(this_a_matrix)
    mr_document["files"]["nv53"]["name"] = f"a{this_cluster_id}_dmtrx.tsv"
    mr_document["files"]["nv53"]["blob"] = this_a_matrix
    mr_document["files"]["bm00"]["name"] = f"b{this_cluster_id}_dmtrx.tsv"
    mr_document["files"]["bm00"]["blob"] = this_a_matrix # TODO: replace with this_b_matrix!!
    mr_document['panes']['model']['layout']['children'][0]['children'][0]['children'][1]['children'][0]['name'] = "Matrix"
    mr_document['panes']['model']['layout']['children'][0]['children'][0]['children'][1]['children'][1]['name'] = "Matrix-copy"

    # actually upload
    assert URL is not None, f"No Microreact URL for {this_cluster_id}!"
    update_resp = requests.post("https://microreact.org/api/projects/update",
        headers={"Access-Token": token, "Content-Type": "application/json; charset=UTF-8"},
        params={"project": URL, "access": "private"},
        timeout=100,
        json=mr_document)
    if update_resp.status_code == 200:
        URL = update_resp.json()['id']
        if verbose: print(f"Updated MR project with id {URL}")
    else:
        print(f"Failed to update MR project with id {URL}:", update_resp.status_code, update_resp.text)

print(all_cluster_information)
os.remove(args.persistentclustermeta)
os.remove(args.persistentids)
all_cluster_information.write_ndjson('all_cluster_information.json')
new_persistent_meta = all_cluster_information.select(['cluster_id', 'first_found', 'last_update', 'jurisdictions', 'microreact_url'])
new_persistent_meta.write_csv(f'persistentMETA{today}.tsv', separator='\t')
new_persistent_ids = hella_redundant.select(['sample_id', 'cluster_id', 'cluster_distance'])
new_persistent_ids.write_csv(f'persistentIDS{today}.tsv', separator='\t')
