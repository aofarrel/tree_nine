VERSION = "0.2.0"
verbose = True
print(f"PROCESS CLUSTERS - VERSION {VERSION}")

# pylint: disable=too-many-statements,too-many-branches,simplifiable-if-expression,too-many-locals,too-complex,consider-using-tuple
# pylint: disable=wrong-import-position,unspecified-encoding,useless-suppression,multiple-statements,line-too-long

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
import logging
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
today = datetime.utcnow().date() # I don't care if this runs past midnight, give everything the same day!
print(f"It's {today} in Thurles right now. Up Tipp!")

def main():
    parser = argparse.ArgumentParser(description="Crunch data, extract trees, upload to MR, etc")
    parser.add_argument('-s', '--shareemail', type=str, required=False, help="email (just one) for calling MR share API")
    parser.add_argument('-to', '--token', type=str, required=True, help="TXT: MR token")
    parser.add_argument('-ls', '--latestsamples', type=str, help='TSV: latest sample information')
    #parser.add_argument('-sm', '--samplemeta', type=str, help='TSV: sample metadata pulled from terra (including myco outs), one line per sample')
    parser.add_argument('-pcm', '--persistentclustermeta', type=str, help='TSV: persistent cluster metadata from last full run of TB-D')
    #parser.add_argument('-ccm', '--latestclustermeta', type=str, help='TSV: latest cluster metadata (as identified by find_clusters.py like five minutes ago)')
    parser.add_argument('-pid', '--persistentids', type=str, help='TSV: persistent IDs from last full run of TB-D')
    #parser.add_argument('-cid', '--latestids', type=str, help='TSV: latest cluster IDs (as identified by find_clusters.py like five minutes ago)')
    #parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')
    parser.add_argument('-mat', '--mat_tree', type=str, help='PB: tree')
    #parser.add_argument('-cs', '--contextsamples', type=int, default=0, help="[UNUSED] int: Number of context samples for cluster subtrees")
    parser.add_argument('-cd', '--combineddiff', type=str, help='diff: Maple-formatted combined diff file, needed for backmasking')
    parser.add_argument('-dl', '--denylist', type=str, required=False, help='TXT: newline delimited list of cluster IDs to never use')

    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO)
    all_latest_samples = pl.read_csv(args.latestsamples, separator="\t", schema_overrides={"latest_cluster_id": pl.Utf8})
    all_persistent_samples = pl.read_csv(args.persistentids, separator="\t", schema_overrides={"cluster_id": pl.Utf8})
    persistent_clusters_meta = pl.read_csv(args.persistentclustermeta, separator="\t", null_values="NULL", try_parse_dates=True, schema_overrides={"cluster_id": pl.Utf8})
    #latest_clusters = pl.read_csv(args.latestclustermeta, separator="\t")
    with open(args.token, 'r') as file:
        try:
            token = file.readline()
        except UnicodeDecodeError:
            token = None # TODO: we're currently only doing this to allow us to run most of the script but not upload, but in prod we'd want to fail ASAP

    print_df_to_debug_log("all_latest_samples", all_latest_samples)
    print_df_to_debug_log("all_persistent_samples", all_persistent_samples)

    # ensure each sample in latest-clusters has, at most, one 20 SNP, one 10 SNP, and one 05 SNP
    temp_latest_groupby_sample = all_latest_samples.group_by("sample_id", maintain_order=True).agg(pl.col("cluster_distance"))
    for row in temp_latest_groupby_sample.iter_rows(named=True):
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
            logging.error(f"{row['sample_id']} has invalid clusters: {row['cluster_distance']}") #pylint: disable=logging-fstring-interpolation
            raise ValueError
    
    # ensure that we didn't have a total swap of cluster names due to samples being removed
    temp_latest_groupby_cluster = all_latest_samples.group_by("latest_cluster_id", maintain_order=True).agg(pl.col("sample_id")).rename({"latest_cluster_id": "cluster_id"})
    temp_persis_groupby_cluster = all_persistent_samples.group_by("cluster_id", maintain_order=True).agg(pl.col("sample_id"))
    print_df_to_debug_log("Latest IDs, grouped by cluster", temp_latest_groupby_cluster)
    print_df_to_debug_log("Persistent IDs, grouped by cluster", temp_persis_groupby_cluster)
    existing_cluster_ids = set(temp_latest_groupby_cluster["cluster_id"].to_list()) | set(temp_persis_groupby_cluster["cluster_id"].to_list())
    latest_overrides = {}
    for latest_row in temp_latest_groupby_cluster.iter_rows(named=True):
        cluster_id = latest_row["cluster_id"]
        latest_samps = set(latest_row["sample_id"])
        persis_row = temp_persis_groupby_cluster.filter(pl.col("cluster_id") == cluster_id)

        if persis_row.is_empty():
            logging.debug("%s present in latest but not persistent", cluster_id)
        else:
            persis_samps = set(persis_row["sample_id"].item())
            if latest_samps.isdisjoint(persis_samps):
                logging.warning("Disjoint union of %s:", cluster_id)
                logging.warning("    This run, %s had these samples: %s", cluster_id, latest_samps)
                logging.warning("    But previous run, %s had these samples: %s", cluster_id, persis_samps)

                # assign a new cluster_id to the latest cluster to avoid reusing a persistent ID
                new_cluster_id = generate_truly_unique_cluster_id(existing_cluster_ids, args.denylist)
                latest_overrides[cluster_id] = str(new_cluster_id).zfill(6)
                existing_cluster_ids.add(new_cluster_id)
                logging.warning("Generated new cluster ID: %s → %s", cluster_id, new_cluster_id)
                subprocess.run(f"mv a{cluster_id}.nwk a{new_cluster_id}.nwk", shell=True, check=True)
                logging.warning("Renamed nwk")
                subprocess.run(f"mv a{cluster_id}.pb a{new_cluster_id}.pb", shell=True, check=True)
                logging.warning("Renamed pb")
                subprocess.run(f"mv a{cluster_id}_dmtrx.tsv a{new_cluster_id}_dmtrx.tsv", shell=True, check=True)
                logging.warning("Renamed distance matrix")
    
    if latest_overrides:
        all_latest_samples = all_latest_samples.with_columns(
            pl.when(pl.col("latest_cluster_id").is_in(list(latest_overrides.keys())))
            .then(pl.col("latest_cluster_id").replace(latest_overrides))
            .otherwise(pl.col("latest_cluster_id"))
            .alias("latest_cluster_id")
        )
        print_df_to_debug_log("Latest IDs after name changes", all_latest_samples)

        # write to denylist so subsequent runs will still avoid this cluster name
        with open('clusterid_denylist.txt', 'w') as file:
            for key in latest_overrides:
                file.write(str(key) + '\n')
    else:
        logging.debug("Did not change any persistent ID names prior to running the main script that also changes persistent IDs (just roll with it)")


    # cluster IDs @ 20, 10, and 5 to prepare for persistent cluster ID assignments
    all_latest_20  = all_latest_samples.filter(pl.col("cluster_distance") == 20).select(["sample_id", "latest_cluster_id"])
    all_latest_10  = all_latest_samples.filter(pl.col("cluster_distance") == 10).select(["sample_id", "latest_cluster_id"])
    all_latest_5   = all_latest_samples.filter(pl.col("cluster_distance") == 5).select(["sample_id", "latest_cluster_id"])
    all_latest_unclustered = all_latest_samples.filter(pl.col("cluster_distance") == -1).select(["sample_id", "latest_cluster_id"]) # pylint: disable=unused-variable

    all_persistent_20 = all_persistent_samples.filter(pl.col("cluster_distance") == 20).select(["sample_id", "cluster_id"])
    all_persistent_10  = all_persistent_samples.filter(pl.col("cluster_distance") == 10).select(["sample_id", "cluster_id"])
    all_persistent_5   = all_persistent_samples.filter(pl.col("cluster_distance") == 5).select(["sample_id", "cluster_id"])
    all_persistent_unclustered = all_persistent_samples.filter(pl.col("cluster_distance") == -1).select(["sample_id", "cluster_id"]) # pylint: disable=unused-variable

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
    print_df_to_debug_log("filtered_latest_20", filtered_latest_20)
    print_df_to_debug_log("filtered_latest_10", filtered_latest_10)
    print_df_to_debug_log("filtered_latest_5", filtered_latest_5)

    filtered_persistent_20 = all_persistent_20.join(all_latest_20.drop(['latest_cluster_id']), on="sample_id", how="inner")
    filtered_persistent_10 = all_persistent_10.join(all_latest_10.drop(['latest_cluster_id']), on="sample_id", how="inner")
    filtered_persistent_5 = all_persistent_5.join(all_latest_5.drop(['latest_cluster_id']), on="sample_id", how="inner")
    print_df_to_debug_log("filtered_persistent_20", filtered_persistent_20)
    print_df_to_debug_log("filtered_persistent_10", filtered_persistent_10)
    print_df_to_debug_log("filtered_persistent_5", filtered_persistent_5)

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

    rosetta_20 = pl.read_csv("rosetta_stone_20.tsv", separator="\t", has_header=False, schema_overrides={"column_2": pl.Utf8}).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})
    rosetta_10 = pl.read_csv("rosetta_stone_10.tsv", separator="\t", has_header=False, schema_overrides={"column_2": pl.Utf8}).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})
    rosetta_5 = pl.read_csv("rosetta_stone_5.tsv", separator="\t", has_header=False, schema_overrides={"column_2": pl.Utf8}).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id'})

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

    print_df_to_debug_log("latest_samples_translated before pl.coalesce", latest_samples_translated)

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

    print_df_to_debug_log("latest_samples_translated after pl.coalesce and check", latest_samples_translated)

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
        logging.info(grouped)
        raise ValueError("Some clusters have multiple unique cluster_distance values.")

    grouped = grouped.with_columns(
        grouped["distance_values"].list.get(0).alias("cluster_distance")
    ).drop(["distance_nunique", "distance_values"])

    print_df_to_debug_log("After grouping and then intager-a-fy", grouped)

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
        # drop AFTER this logging.info so we can see if anything is whacky here
        logging.info("After looking for new samples")
        logging.info(grouped)

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
    print_df_to_debug_log("Joined grouped with latest_samples_translated to form hella_redundant, and sorted by cluster distance", hella_redundant)

    # this is a super goofy way to link parents and children, but it seems to work
    # basically, we're building a dataframe that has per-sample information, and a bunch of
    # duplicated per-cluster information. this is EXTREMELY redundant and in a normal world
    # we would not do this!
    sample_map = {dist: {} for dist in [5, 10, 20]}
    for row in hella_redundant.iter_rows(named=True):
        sample_map[row["cluster_distance"]][row["sample_id"]] = row["cluster_id"]
    logging.debug("Sample map:")
    logging.debug(sample_map)
    updates = []
    for row in hella_redundant.iter_rows(named=True):
        cluster_id, one_sample, distance = row["cluster_id"], row["sample_id"], row["cluster_distance"]
        logging.info("[{%i}] sample %s in cluster_id %s", distance, one_sample, cluster_id)
        if distance == 5:
            parent_id = sample_map[10].get(one_sample)
            if parent_id:
                updates.append((cluster_id, "cluster_parent", parent_id))
        elif distance == 10:
            parent_id = sample_map[20].get(one_sample)
            if parent_id:
                updates.append((cluster_id, "cluster_parent", parent_id))
            child_id = sample_map[5].get(one_sample)  # INTENTIONALLY ONLY GRABS THE CHILD ID ASSOCIATED WITH THIS SAMPLE, NOT ALL CHILDREN OF CLUSTER
            if child_id:
                updates.append((cluster_id, "cluster_one_child", child_id))
        elif distance == 20:
            child_id = sample_map[10].get(one_sample)  # INTENTIONALLY ONLY GRABS THE CHILD ID ASSOCIATED WITH THIS SAMPLE, NOT ALL CHILDREN OF CLUSTER
            if child_id:
                updates.append((cluster_id, "cluster_one_child", child_id))
        else:
            raise ValueError
    print_df_to_debug_log("updates", updates)
    for cluster_id, col, value in updates:
        logging.debug("For cluster %s, col %s, val %s in updates", cluster_id, col, value)
        if col == "cluster_parent":
            hella_redundant = update_cluster_column(hella_redundant, cluster_id, "cluster_parent", value)
        else:
            hella_redundant = hella_redundant.with_columns(
                pl.when(pl.col("cluster_id") == cluster_id)
                .then((pl.col("cluster_children").list.concat(pl.lit(value))).list.unique())
                .otherwise(pl.col("cluster_children"))
                .alias("cluster_children")
            )
    cluster_id = None
    print_df_to_debug_log("after linking parents and children", hella_redundant)

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

    print_df_to_debug_log("after processing what clusters are brand new and what ones are unchanged", hella_redundant)

    # drop AFTER the verbose logging.info
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
    print_df_to_debug_log("after grouping hella_redundant by cluster_id", second_group)

    # check cluster distances
    # TODO: this and other asserts will probably need to change if we change how we handle unclustered samples
    assert (second_group["cluster_distance"].list.len() == 1).all(), "cluster_distance lists have length ≠ 1"
    second_group = second_group.with_columns(pl.col("cluster_distance").list.get(0).alias("cluster_distance_int"))
    second_group = second_group.drop("cluster_distance").rename({"cluster_distance_int": "cluster_distance"})
    logging.info(second_group)

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
    all_cluster_information = get_nwk_and_matrix_plus_local_mask(all_cluster_information, args.combineddiff)

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
        all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("last_update"))
    if "first_found" not in all_cluster_information.columns:
        all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("first_found"))

    print_df_to_debug_log("Old persistent cluster metadata", persistent_clusters_meta)
    print_df_to_debug_log("All cluster information we have so far after joining with that", all_cluster_information)

    # okay, everything looks good so far. let's get some URLs!!
    logging.info("Assigning self-URLs...")
    for row in all_cluster_information.iter_rows(named=True):
        this_cluster_id = row["cluster_id"]
        distance = row["cluster_distance"]
        has_children = False if row["cluster_children"] is None else True
        brand_new = row["cluster_brand_new"]
        needs_updating = row["cluster_needs_updating"]
        URL = row["microreact_url"]

        # if a childless 20-cluster is brand new, it gets a found and an update date, but no MR URL
        # if a childless 20-cluster isn't new but has new samples, we change the update date (but still no MR URL)

        if brand_new:
            assert URL is None, f"{this_cluster_id} is brand new but already has a MR URL?"
            all_cluster_information = update_first_found(all_cluster_information, this_cluster_id)
            all_cluster_information = update_last_update(all_cluster_information, this_cluster_id)
            if distance == 20 and not has_children:
                logging.info("%s is a brand-new 20-cluster with no children, not uploading", this_cluster_id)
                all_cluster_information = update_first_found(all_cluster_information, this_cluster_id)
                all_cluster_information = update_last_update(all_cluster_information, this_cluster_id)
                all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "cluster_needs_updating", False)
                continue
            URL = create_new_mr_project(token, this_cluster_id)
            all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "microreact_url", URL)

        elif needs_updating:
            # later on, assert URL and first_found is not None... but for the first time don't do that!
            if URL is None:
                if distance == 20 and not has_children:
                    logging.info("%s is an old 20-cluster with no children, not uploading", this_cluster_id)
                    all_cluster_information = update_last_update(all_cluster_information, this_cluster_id)
                    all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "cluster_needs_updating", False)
                    continue
                logging.warning("%s isn't brand new, but is flagged as needing an update and has no URL. Will make a new URL.", this_cluster_id)
                URL = create_new_mr_project(token, this_cluster_id)
                all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "microreact_url", URL)
                all_cluster_information = update_last_update(all_cluster_information, this_cluster_id)
            else:
                logging.info("%s's URL (%s) seems valid, but we're not gonna check it because that'd be a pain in the neck", this_cluster_id, URL)

        else:
            # unchanged
            continue

    all_cluster_information = all_cluster_information.with_columns(
        pl.lit(None).alias("parent_URL"),
        pl.lit(None).alias("children_URLs") # TODO: how are gonna keep these ordered with the children column... does that even matter?
    )

    # now that everything has a URL, or doesn't need one, iterate a second time to get URLs of parents and children
    logging.info("Searching for MR URLs of parents and children...")
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
            if row["sample_id"] is not None:
                logging.debug("%s is a %i-cluster with no new samples (ergo no new children), skipping", this_cluster_id, distance)
            else:
                logging.warning("%s has no samples! This is likely a decimated cluster that lost all of its samples. We will not be updating its MR project.") # TODO but we should
                continue

        if has_parent:
            # get value of column "microreact_url" when column "cluster_id" matches the string in row["cluster_parent"]
            parent_URL = all_cluster_information.filter(
                all_cluster_information["cluster_id"] == row["cluster_parent"]
            )["microreact_url"].item()
            logging.debug("%s has parent URL %s", this_cluster_id, parent_URL)
            all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "parent_URL", parent_URL)
            
        if has_children:
            # get value of column "microreact_url" when column "cluster_id" matches a string in the list in row["cluster_children"]
            children_URLs = all_cluster_information.filter(
                all_cluster_information["cluster_id"].is_in(row["cluster_children"])
            )["microreact_url"].to_list()
            logging.debug("%s has children with URLs %s", this_cluster_id, children_URLs)
            all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "children_URLs", children_URLs)
    parent_URL, children_URLs, URL = None, None, None
    logging.info("Final(ish) data table:")
    logging.info(all_cluster_information)
    all_cluster_information.write_ndjson('all_cluster_information_before_upload.json')

    # now that everything can be crosslinked, iterate one more time to actually upload
    # yes three iterations is weird, cringe even. I'm doing this to make debugging easier.
    logging.info("Splitting and uploading...")
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
        URL = row["microreact_url"]
        try:
            first_found = today if row["first_found"] is None else datetime.strptime(row["first_found"], "%Y-%m-%d")
            first_found_shorthand = f'{str(first_found.year)}M{str(first_found.month).zfill(2)}'
        except TypeError: # fires if !today and first_found is datetime.date (as opposed to str)
            first_found = today if row["first_found"] is None else str(row["first_found"].isoformat())
            first_found_shorthand = f'{str(row["first_found"].year)}M{str(row["first_found"].month).zfill(2)}'

        # Because there is never a situation where a new child cluster pops up in a parent cluster that doesn't need to be updated,
        # and because MR URLs don't need to be updated, clusters that don't need updating don't need to know parent/child URLs.
        if not needs_updating:
            if row["sample_id"] is not None:
                logging.info("%s is a %i-cluster marked as not needing updating (no new samples and/or childless 20-cluster), skipping", this_cluster_id, distance)
            else:
                logging.warning("%s has no samples! This is likely a decimated cluster that lost all of its samples. We will not be updating its MR project.") # TODO but we should
                continue
        
        with open("./REALER_template.json", "r") as real_template_json:
            mr_document = json.load(real_template_json)

        # project title
        fullID = f"{str(distance).zfill(2)}SNP-{first_found_shorthand}-{this_cluster_id}"
        mr_document["meta"]["name"] = f"{fullID} updated {today.isoformat()}"
        mr_document["meta"]["description"] = f"{this_cluster_id} as automatically generated by version {VERSION} of process_clusters.py"

        # note
        markdown_note = f"### {this_cluster_id} ({distance}-SNP, {len(sample_id_list)} samples)\n*Updated {today.isoformat()}*\n\n"
        #if len(sample_id_list) == 2:
        #    markdown_note += "*WARNING: If this cluster's SNP distances are all 0, it may not render correctly in Microreact*\n\n"
        markdown_note += f"First found {first_found.isoformat()}, UUID {this_cluster_id}, fullID {fullID}\n\n"
        if has_parent:
            markdown_note += f"Parent cluster: [{cluster_parent}](https://microreact.org/project/{parent_URL})\n\n"
        if has_children:
            markdown_note += f"Child clusters ({n_children} total):\n" + "".join(f"* [click here](https://microreact.org/project/{child_url})\n" for child_url in row["children_URLs"]) + "\n"
        markdown_note += f"\n\nCluster-first-found and last-update dates are calculated as UST when the pipeline was run; dates in the metadata table are untouched. process_clusters.py version: {VERSION}\n"
        mr_document["notes"]["note-1"]["source"] = markdown_note

        # trees
        this_a_nwk = get_atree_raw(this_cluster_id, all_cluster_information)
        this_b_nwk = get_btree_raw(this_cluster_id, all_cluster_information)
        mr_document["files"]["chya"]["name"] = f"a{this_cluster_id}.nwk" # can be any arbitrary name since we already extracted nwk
        mr_document["files"]["chya"]["blob"] = this_a_nwk
        mr_document["files"]["bmtr"]["name"] = f"b{this_cluster_id}.nwk"
        mr_document["files"]["bmtr"]["blob"] = this_b_nwk
        mr_document['panes']['model']['layout']['children'][0]['children'][0]['children'][0]['children'][0]['name'] = "Raw Tree"
        mr_document['panes']['model']['layout']['children'][0]['children'][0]['children'][0]['children'][1]['name'] = "Locally Masked (Bionumerics-style)"

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
        logging.debug(labels)
        mr_document["files"]["ji0o"]["name"] = f"{this_cluster_id}_metadata.csv"
        mr_document["files"]["ji0o"]["blob"] = labels

        # distance matrix already readlines()'d into this_a_matrix, just make it a string
        this_a_matrix = get_amatrix_raw(this_cluster_id, all_cluster_information)
        this_a_matrix = "\n".join(this_a_matrix)
        this_b_matrix = get_bmatrix_raw(this_cluster_id, all_cluster_information)
        this_b_matrix = "\n".join(this_b_matrix)
        mr_document["files"]["nv53"]["name"] = f"a{this_cluster_id}_dmtrx.tsv"
        mr_document["files"]["nv53"]["blob"] = this_a_matrix
        mr_document["files"]["bm00"]["name"] = f"b{this_cluster_id}_dmtrx.tsv"
        mr_document["files"]["bm00"]["blob"] = this_b_matrix
        mr_document['panes']['model']['layout']['children'][0]['children'][0]['children'][1]['children'][0]['name'] = "Raw Matrix"
        mr_document['panes']['model']['layout']['children'][0]['children'][0]['children'][1]['children'][1]['name'] = "Locally Masked (Bionumerics-style)"

        # actually upload
        assert URL is not None, f"No Microreact URL for {this_cluster_id}!"
        update_resp = requests.post("https://microreact.org/api/projects/update",
            headers={"Access-Token": token, "Content-Type": "application/json; charset=UTF-8"},
            params={"project": URL, "access": "private"},
            timeout=100,
            json=mr_document)
        if update_resp.status_code == 200:
            URL = update_resp.json()['id']
            logging.debug("Updated MR project with id %s", URL)

            # share the project
            #share_mr_project(token, URL, emails) 

        else:
            logging.error("Failed to update MR project with id %s [code %s]: %s", URL, update_resp.status_code, update_resp.text) 
            logging.error("Will continue...")

    logging.info("Finished. All cluster information:")
    logging.info(all_cluster_information)
    os.remove(args.persistentclustermeta)
    os.remove(args.persistentids)
    os.remove(args.token)
    all_cluster_information.write_ndjson('all_cluster_information.json')
    new_persistent_meta = all_cluster_information.select(['cluster_id', 'first_found', 'last_update', 'jurisdictions', 'microreact_url'])
    new_persistent_meta.write_csv(f'persistentMETA{today.isoformat()}.tsv', separator='\t')
    new_persistent_ids = hella_redundant.select(['sample_id', 'cluster_id', 'cluster_distance'])
    new_persistent_ids.write_csv(f'persistentIDS{today.isoformat()}.tsv', separator='\t')


def add_col_if_not_there(dataframe, column):
    if column not in dataframe.columns:
        return dataframe.with_columns(pl.lit(None).alias(column))
    return dataframe

def get_nwk_and_matrix_plus_local_mask(big_ol_dataframe, combineddiff):
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "a_matrix")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "a_tree")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_matrix")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_tree")
    print("get_nwk_and_matrix_plus_local_mask() got this dataframe", big_ol_dataframe)
    for row in big_ol_dataframe.iter_rows(named=True):
        this_cluster_id = row["cluster_id"]
        workdir_cluster_id = row["workdir_cluster_id"]
        amatrix = f"a{workdir_cluster_id}_dmtrx.tsv" if os.path.exists(f"a{workdir_cluster_id}_dmtrx.tsv") else None
        atree = f"a{workdir_cluster_id}.nwk" if os.path.exists(f"a{workdir_cluster_id}.nwk") else None

        logging.debug("[%s] found %s and %s", this_cluster_id, amatrix, atree)
        big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "a_matrix", amatrix)
        big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "a_tree", atree)
        logging.debug("[%s] updated df, now tryna rename files if necessary", this_cluster_id)
        if workdir_cluster_id != this_cluster_id: # do NOT remove this check
            if amatrix is not None:
                if not os.path.exists(f"a{this_cluster_id}_dmtrx.tsv"): # and that's why we can't remove aforementioned check
                    os.rename(f"a{workdir_cluster_id}_dmtrx.tsv", f"a{this_cluster_id}_dmtrx.tsv")
                    amatrix = f"a{this_cluster_id}_dmtrx.tsv"
                else:
                    logging.debug("[%s] Cannot rename a{%s}_dmtrx.tsv to a{%s}_dmtrx.tsv, will rename to a{%s}_dmtrx_temp.tsv for now", this_cluster_id, workdir_cluster_id, this_cluster_id, this_cluster_id)
                    os.rename(f"a{workdir_cluster_id}_dmtrx.tsv", f"a{this_cluster_id}_dmtrx_temp.tsv")
                    amatrix = f"a{this_cluster_id}_dmtrx_temp.tsv" # TODO: after done iterating, we should be able to remove _temp from files and update df accordingly
                big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "a_matrix", amatrix)

            if atree is not None:
                if not os.path.exists(f"a{this_cluster_id}.nwk"):
                    os.rename(f"a{workdir_cluster_id}.nwk", f"a{this_cluster_id}.nwk")
                    atree = f"a{this_cluster_id}.nwk"
                else:
                    logging.debug("[%s] Cannot rename a{%s}.nwk to a{%s}.nwk, will rename to a{%s}_temp.nwk for now", this_cluster_id, workdir_cluster_id, this_cluster_id, this_cluster_id)
                    os.rename(f"a{workdir_cluster_id}.nwk", f"a{this_cluster_id}_temp.nwk")
                    atree = f"a{this_cluster_id}_temp.nwk" # TODO: after done iterating, we should be able to remove _temp from files and update df accordingly
                big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "a_tree", atree)
        
        logging.debug("[%s] now dealing with the b-sides", this_cluster_id)
        btree = bmatrix = None
        if atree is not None:
            logging.debug("[%s] atree is not none", this_cluster_id)
            atreepb = next((f"a{id}.pb" for id in [this_cluster_id, workdir_cluster_id] if os.path.exists(f"a{id}.pb")), None)
            if atreepb:
                logging.debug("[%s] atreepb is not none", this_cluster_id)
                btreepb = f"b{this_cluster_id}.pb"
                btree = f"b{this_cluster_id}.nwk"
                try:
                    subprocess.run(f"matUtils mask -i {atreepb} -o {btreepb} -D 1000 -f {combineddiff}", shell=True, check=True)
                    logging.debug("[%s] matUtils mask returned 0 (atree.pb --> masked btree.pb)", this_cluster_id)
                    subprocess.run(f"matUtils extract -i {btreepb} -t {btree}", shell=True, check=True)
                    logging.debug("[%s] matUtils extract returned 0 (masked btree.pb --> masked btree.nwk)", this_cluster_id)
                    subprocess.run(f"python3 /scripts/find_clusters.py {btreepb} {btree} --type BM --collection-name {this_cluster_id} --nocluster --nolonely --noextraouts", shell=True, check=False)
                    logging.debug("[%s] ran find_clusters.py but DID NOT CHECK EXIT CODE so there's a chance bmatrix doesn't exist", this_cluster_id)
                    bmatrix = f"b{this_cluster_id}_dmtrx.tsv" if os.path.exists(f"b{this_cluster_id}_dmtrx.tsv") else None
                except subprocess.CalledProcessError as e:
                    logging.warning("[%s] Failed to generate locally-masked tree/matrix: %s", this_cluster_id, e.output)
                big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_matrix", bmatrix)
                big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_tree", btree)
    print_df_to_debug_log("returning this dataframe from get_nwk_and_matrix_plus_local_mask()", big_ol_dataframe)
    return big_ol_dataframe

def update_cluster_column(df, cluster_id, column, new_value):
    assert column in df.columns, f"Tried to update {column} with {new_value} but {column} not in dataframe?"
    return df.with_columns(
        pl.when(df["cluster_id"] == cluster_id)
        .then(pl.lit(new_value))
        .otherwise(df[column])
        .alias(column))

def update_first_found(df, cluster_id):
    return df.with_columns(
        pl.when(df["cluster_id"] == cluster_id)
        .then(pl.lit(today.isoformat()))
        .otherwise(df["first_found"])
        .alias("first_found"))

def update_last_update(df, cluster_id):
    return df.with_columns(
        pl.when(df["cluster_id"] == cluster_id)
        .then(pl.lit(today.isoformat()))
        .otherwise(df["last_update"])
        .alias("last_update"))

def print_df_to_debug_log(dataframe_name, actual_dataframe):
    logging.debug("%s", dataframe_name)
    logging.debug(actual_dataframe)

def share_mr_project(token, mr_url, email):
    api_url = "https://microreact.org/api/shares/add"
    params = {"id": mr_url}
    headers = {"Access-Token": token}
    data = {"emails": [email], "role": "viewer" }
    share_resp = requests.post(api_url, headers=headers, params=params, json=data, timeout=100)
    if share_resp.status_code == 200: 
        logging.debug("Successfully shared project")
    else:
        logging.error("Failed to share MR project with id %s [code %s]: %s", mr_url, share_resp.status_code, share_resp.text)

def create_new_mr_project(token, this_cluster_id):
    with open("./BLANK_template.json", "r") as temp_proj_json:
        mr_document = json.load(temp_proj_json)
    update_resp = requests.post("https://microreact.org/api/projects/create",
        headers={"Access-Token": token, "Content-Type": "application/json; charset=UTF-8"},
        params={"access": "private"},
        timeout=100,
        json=mr_document)
    if update_resp.status_code == 200:
        URL = update_resp.json()['id']
        logging.debug("Created temporary MR project with id %s", URL)
    else:
        logging.warning("Failed to create new MR project [%s]: %s", update_resp.status_code, update_resp.text)
    logging.info("%s is brand new and has been assigned MR URL %s and a first-found date", this_cluster_id, URL)
    return URL

def get_atree_raw(cluster_name, big_ol_dataframe):
    try:
        atree_series = big_ol_dataframe.filter(pl.col("cluster_id") == cluster_name).select("a_tree")
        atree = atree_series.item()
        with open(atree, "r") as nwk_file:
            return nwk_file.readline() # only need first line
    except (OSError, TypeError): # OSError: File Not Found, TypeError: None
        return "((INITIAL_SUBTREE_ERROR:1,REPORT_THIS_BUG_TO_ASH:1):1,DO_NOT_INCLUDE_PHI_IN_REPORT:1);"

def get_btree_raw(cluster_name, big_ol_dataframe):
    try:
        btree_series = big_ol_dataframe.filter(pl.col("cluster_id") == cluster_name).select("b_tree")
        btree = btree_series.item()
        with open(btree, "r") as nwk_file:
            return nwk_file.readline() # only need first line
    except (OSError, TypeError):
        return "((MASKED_SUBTREE_ERROR:1,REPORT_THIS_BUG_TO_ASH:1):1,DO_NOT_INCLUDE_PHI_IN_REPORT:1);"

def generate_truly_unique_cluster_id(existing_ids, denylist):
    if denylist is not None:
        with open(denylist, "r") as f:
            denylist = {line.strip() for line in f}
    else:
        denylist = set()
    existing_set = set(existing_ids) | set(denylist)
    new_id = max((int(x) for x in existing_ids if x.isdigit()), default=0) + 1  # Start from max valid int + 1
    while str(new_id).zfill(6) in existing_set:
        new_id += 1
    return str(new_id).zfill(6)

def get_amatrix_raw(cluster_name, big_ol_dataframe):
    try:
        amatrix_series = big_ol_dataframe.filter(pl.col("cluster_id") == cluster_name).select("a_matrix")
        amatrix = amatrix_series.item()
        with open(amatrix, "r") as distance_matrix:
            this_a_matrix = distance_matrix.readlines()
        return this_a_matrix
    except (OSError, TypeError):
        return """INITIAL_SUBTREE_ERROR\tREPORT_THIS_BUG_TO_ASH\tDO_NOT_INCLUDE_PHI_IN_REPORT\n
        INITIAL_SUBTREE_ERROR\t1\t1\t1\n
        REPORT_THIS_BUG_TO_ASH\t1\t1\t1\n
        DO_NOT_INCLUDE_PHI_IN_REPORT\t1\t1\t1"""

def get_bmatrix_raw(cluster_name, big_ol_dataframe):
    try:
        bmatrix_series = big_ol_dataframe.filter(pl.col("cluster_id") == cluster_name).select("b_matrix")
        bmatrix = bmatrix_series.item()
        with open(bmatrix, "r") as distance_matrix:
            this_b_matrix = distance_matrix.readlines()
        return this_b_matrix
    except (OSError, TypeError):
        return """MASKED_SUBTREE_ERROR\tREPORT_THIS_BUG_TO_ASH\tDO_NOT_INCLUDE_PHI_IN_REPORT\n
        MASKED_SUBTREE_ERROR\t1\t1\t1\n
        REPORT_THIS_BUG_TO_ASH\t1\t1\t1\n
        DO_NOT_INCLUDE_PHI_IN_REPORT\t1\t1\t1"""

if __name__ == "__main__":
    main()