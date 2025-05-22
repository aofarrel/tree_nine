VERSION = "0.3.3" # does not necessarily match Tree Nine git version
verbose = True
print(f"PROCESS CLUSTERS - VERSION {VERSION}")

# pylint: disable=too-many-statements,too-many-branches,simplifiable-if-expression,too-many-locals,too-complex,consider-using-tuple
# pylint: disable=wrong-import-position,unspecified-encoding,useless-suppression,multiple-statements,line-too-long,consider-using-sys-exit

# Note to future maintainers: We are using "polars" for dataframes here, which is like pandas, but significantly more efficient.
# Based on my experience working with Literally Every Mycobacterium Sample On NCBI SRA's And Its Metadata, I estimate this script
# will be performant on reasonable hardware (eg, Intel-hardware 2019 16 GB MacBook Pro) up until *roughly* 50,000 samples.
#
# This is not to say this script is optimized. It creates several extremely redundant dataframes. You will see dataframes
# be joined, aggregated, then exploded, all over the place. We also use sample-level dataframes several times, where
# each sample is present up to three times, because every sample can be a member of a 20 SNP cluster, a 10 SNP cluster,
# and a 5 SNP cluster. While there are some parts that could likely be made more effecient, the overall "structure" of
# handling 20/10/5 SNPs seperately until we're ready to make parent-child connections is probably for the best, as is
# using sample-level dataframes as necessary, redundant as they may be. In other words -- I do not recommend attempting
# a major refactor of this script, nor do I recommend replacing Marc's script with your own solution, because there are
# several edge cases to account for, and because you're unlikely to make meaningful performance gains.
# 

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
from polars.testing import assert_series_equal
pl.Config.set_tbl_rows(-1)
pl.Config.set_tbl_cols(-1)
pl.Config.set_tbl_width_chars(2000)
pl.Config.set_fmt_str_lengths(500)
pl.Config.set_fmt_table_cell_list_len(500)
today = datetime.utcnow().date() # I don't care if this runs past midnight, give everything the same day!
print(f"It's {today} in Thurles right now. Up Tipp!")
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("requests").setLevel(logging.WARNING)

def main():
    parser = argparse.ArgumentParser(description="Crunch data, extract trees, upload to MR, etc")
    parser.add_argument('-s', '--shareemail', type=str, required=False, help="email (just one) for calling MR share API")
    parser.add_argument('-to', '--token', type=str, required=False, help="TXT: MR token")
    parser.add_argument('-as', '--allsamples', type=str, required=False, help='')
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
    parser.add_argument('-mr', '--yes_microreact', action='store_true')
    parser.add_argument('-d', '--today', type=str, required=True, help='ISO 8601 date, YYYY-MM-DD')

    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO)
    all_latest_samples = pl.read_csv(args.latestsamples,
        separator="\t", 
        schema_overrides={"latest_cluster_id": pl.Utf8}).filter(pl.col("latest_cluster_id").is_not_null())
    all_persistent_samples = pl.read_csv(args.persistentids,
        separator="\t", 
        schema_overrides={"cluster_id": pl.Utf8}).filter(pl.col("cluster_id").is_not_null())
    persistent_clusters_meta = pl.read_csv(args.persistentclustermeta,
        separator="\t",
        null_values="NULL",
        try_parse_dates=True, 
        schema_overrides={"cluster_id": pl.Utf8}).filter(pl.col("cluster_id").is_not_null())
    #latest_clusters_meta = pl.read_csv(args.latestclustermeta, separator="\t")

    global today # pylint: disable=global-statement
    args_today = datetime.strptime(args.today, "%Y-%m-%d").date()
    if args_today != today:
        # this is just to warn the user they might be using an old or cached date, but we have
        # to use the user-provided date anyway for WDL output matching to work without glob()
        # (we are avoiding WDL-glob() because it creates a random-name folder in GCP which is annoying)
        logging.warning("The date you provided (%s, interpreted as type %s) doesn't match the date in Thurles.", args_today, type(args_today))
        today = args_today

    if args.yes_microreact and not args.token:
        logging.error("You entered --yes_microreact but didn't provide a token file with --token")
        raise ValueError
    if args.token:
        with open(args.token, 'r') as file:
            token = file.readline()
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
    
    persis_groupby_cluster = all_persistent_samples.group_by("cluster_id", maintain_order=True).agg(pl.col("sample_id"), pl.col("cluster_distance").unique().first())
    # This is very tricky: We need to make sure that if any persistent clusters don't exist anymore, their IDs do not get reused.
    # It should only happen when running on a subset of samples and/or if samples have been removed. (In theory something like
    # this can also happen if two clusters get merged, but we'll cross that bridge later.)
    # Previously we would iterate temp_latest_groupby_cluster rowwise and check if
    # this_clusters_latest_samps.isdisjoint(persis_groupby_cluster.filter(pl.col("cluster_id") == latest_cluster_id)),
    # but this was a terrible idea because it would fire whenever samples simply generated a different cluster ID.
    # We instead want to start with a simple question:
    # Are there any samples present in all_persistent_samples not present in all_latest_samples?
    # If no: Literally who cares, the perl script will handle it
    # If yes: Iterate the *persistent* clusters rowwise to make sure they aren't decimated
    all_latest_samples_set = set(all_latest_samples["sample_id"].to_list())
    all_persistent_samples_set = set(all_persistent_samples["sample_id"].to_list())
    logging.debug("Set of all latest samples")
    logging.debug(all_latest_samples_set)
    logging.debug("Set of all persistent samples")
    logging.debug(all_persistent_samples_set)
    if all_persistent_samples_set.issubset(all_latest_samples_set):
        logging.info("All persistent samples is a subset of all latest samples")
    else:
        samples_missing_from_latest = all_persistent_samples_set - all_latest_samples_set # these are sets so this excludes samples exclusive to all_latest
        logging.warning("Samples appear to be missing from the latest run: %s", samples_missing_from_latest)
        if args.allsamples:
            all_input_samples_including_unclustered = args.allsamples.split(',')
        else:
            all_input_samples_including_unclustered = None
            logging.warning("Missing args.allsamples; can't be sure if missing samples are dropped because they no longer cluster or if they were never input.")
        for sample in samples_missing_from_latest:
            if all_input_samples_including_unclustered is None:
                pass
            elif sample in all_input_samples_including_unclustered:
                logging.warning("%s is newly unclustered", sample)
            else:
                logging.warning("%s seems to have been dropped from inputs", sample)
            # get persistent cluster ID regardless
            cluster_ids = get_cluster_ids_for_sample(all_persistent_samples, sample)
            for cluster in cluster_ids:
                if len(get_other_samples_in_cluster(all_persistent_samples, cluster, samples_missing_from_latest)) <= 1:
                    # In theory we could handle this, in practice it's a massive pain in the neck and very easy to mess up!!
                    logging.warning("%s is decimated thanks to losing all samples (or all but one)", cluster)
                    # IF AND ONLY IF this is not on MR (which should only happen if this is a 20-cluster with no subclusters),
                    # we can live with this being a decimated cluster.
                    if not has_microreact_url(persistent_clusters_meta, cluster):
                        logging.info("%s already lacks a Microreact URL, so we can live with it being decimated")
                    else:
                        logging.error("%s has an MR URL and should never be decimated. Cannot continue.")
                        exit(55)
                else:
                    logging.info("Dropped %s from %s but that seems to be okay", sample, cluster)
    
    # cluster IDs @ 20, 10, and 5 to prepare for persistent cluster ID assignments
    all_latest_20  = all_latest_samples.filter(pl.col("cluster_distance") == 20).select(["sample_id", "latest_cluster_id"])
    all_latest_10  = all_latest_samples.filter(pl.col("cluster_distance") == 10).select(["sample_id", "latest_cluster_id"])
    all_latest_5   = all_latest_samples.filter(pl.col("cluster_distance") == 5).select(["sample_id", "latest_cluster_id"])
    all_latest_unclustered = all_latest_samples.filter(pl.col("cluster_distance") == -1).select(["sample_id", "latest_cluster_id"]) # pylint: disable=unused-variable

    all_persistent_20  = all_persistent_samples.filter(pl.col("cluster_distance") == 20).select(["sample_id", "cluster_id"])
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
    filtered_latest_20 = all_latest_20.join(all_persistent_20.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
    filtered_latest_10 = all_latest_10.join(all_persistent_10.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
    filtered_latest_5 = all_latest_5.join(all_persistent_5.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
    print_df_to_debug_log("filtered_latest_20", filtered_latest_20)
    print_df_to_debug_log("filtered_latest_10", filtered_latest_10)
    print_df_to_debug_log("filtered_latest_5", filtered_latest_5)

    filtered_persistent_20 = all_persistent_20.join(all_latest_20.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
    filtered_persistent_10 = all_persistent_10.join(all_latest_10.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
    filtered_persistent_5 = all_persistent_5.join(all_latest_5.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
    print_df_to_debug_log("filtered_persistent_20", filtered_persistent_20)
    print_df_to_debug_log("filtered_persistent_10", filtered_persistent_10)
    print_df_to_debug_log("filtered_persistent_5", filtered_persistent_5)

    filtered_latest_20.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_20.tsv', separator='\t', include_header=False)
    filtered_latest_10.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_10.tsv', separator='\t', include_header=False)
    filtered_latest_5.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_5.tsv', separator='\t', include_header=False)
    filtered_persistent_20.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_20.tsv', separator='\t', include_header=False)
    filtered_persistent_10.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_10.tsv', separator='\t', include_header=False)
    filtered_persistent_5.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_5.tsv', separator='\t', include_header=False)

    subprocess.run("perl /scripts/marcs_incredible_script_update.pl filtered_persistent_20.tsv filtered_latest_20.tsv", shell=True, check=True, capture_output=True, text=True)
    subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_20.tsv", shell=True, check=True)
    subprocess.run("perl /scripts/marcs_incredible_script_update.pl filtered_persistent_10.tsv filtered_latest_10.tsv", shell=True, check=True, capture_output=True, text=True)
    subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_10.tsv", shell=True, check=True)
    subprocess.run("perl /scripts/marcs_incredible_script_update.pl filtered_persistent_5.tsv filtered_latest_5.tsv", shell=True, check=True, capture_output=True, text=True)
    subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_5.tsv", shell=True, check=True)
    if logging.root.level == logging.DEBUG:
        for rock in ['rosetta_stone_20.tsv', 'rosetta_stone_10.tsv', 'rosetta_stone_5.tsv']:
            with open(rock, 'r') as file:
                logging.debug("Contents of %s", rock)
                logging.debug(list(file))
    subprocess.run("/bin/bash /scripts/strip_tsv.sh rosetta_stone_20.tsv rosetta_stone_20_merges.tsv", shell=True, check=True)
    subprocess.run("/bin/bash /scripts/strip_tsv.sh rosetta_stone_10.tsv rosetta_stone_10_merges.tsv", shell=True, check=True)
    subprocess.run("/bin/bash /scripts/strip_tsv.sh rosetta_stone_5.tsv rosetta_stone_5_merges.tsv", shell=True, check=True)
    if logging.root.level == logging.DEBUG:
        for rock in ['rosetta_stone_20.tsv', 'rosetta_stone_10.tsv', 'rosetta_stone_5.tsv', 'rosetta_stone_20_merges.tsv', 'rosetta_stone_10_merges.tsv', 'rosetta_stone_5_merges.tsv']:
            try:
                with open(rock, 'r') as file:
                    subprocess.run(f"/bin/bash /scripts/equalize_tabs.sh {rock}", shell=True, check=True)
                    logging.debug("Contents of %s after splitting and tab-equalizing", rock)
                    logging.debug(list(file))
            except FileNotFoundError:
                logging.debug("Could not find %s but that's probably okay", rock) # can happen if there is no merges

    # we need schema_overrides or else cluster IDs can become non-zfilled i64
    # For some godforesaken reason, some versions of polars will throw `polars.exceptions.ComputeError: found more fields than defined in 'Schema'` even if we set
    # infer_schema = True with a hella large infer_schema_length. Idk why because the exact same file works perfectly fine on my local installation of polars (polars==1.27.0)
    # without even needing to set anything with infer_schema!! Not even a try-except with the except having a three column schema works!! Ugh!!!
    # TODO: is this because the docker is polars==1.26.0?
    # ---> WORKAROUND: equalize_tabs.sh
    rosetta_20 = pl.read_csv("rosetta_stone_20.tsv", separator="\t", has_header=False,
        schema_overrides={"column_1": pl.Utf8, "column_2": pl.Utf8, "column_3": pl.Utf8}, truncate_ragged_lines=True, ignore_errors=True, infer_schema_length=5000).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id', 'column_3': 'special_handling'})
    rosetta_10 = pl.read_csv("rosetta_stone_10.tsv", separator="\t", has_header=False,
        schema_overrides={"column_1": pl.Utf8, "column_2": pl.Utf8, "column_3": pl.Utf8}, truncate_ragged_lines=True, ignore_errors=True, infer_schema_length=5000).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id', 'column_3': 'special_handling'})
    rosetta_5 = pl.read_csv("rosetta_stone_5.tsv", separator="\t", has_header=False, 
        schema_overrides={"column_1": pl.Utf8, "column_2": pl.Utf8, "column_3": pl.Utf8}, truncate_ragged_lines=True, ignore_errors=True, infer_schema_length=5000).rename({'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id', 'column_3': 'special_handling'})

    latest_samples_translated = (all_latest_samples.join(rosetta_20, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_20_cluster_id'}).drop("latest_cluster_id_right")
    latest_samples_translated = (latest_samples_translated.join(rosetta_10, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_10_cluster_id'}).drop("latest_cluster_id_right")
    latest_samples_translated = nullfill_LR(latest_samples_translated, "special_handling", "special_handling_right")
    latest_samples_translated = (latest_samples_translated.join(rosetta_5, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_5_cluster_id'}).drop("latest_cluster_id_right")
    latest_samples_translated = nullfill_LR(latest_samples_translated, "special_handling", "special_handling_right")
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

    print_df_to_debug_log("latest_samples_translated after pl.coalesce and check (sorted by workdir_cluster_id in this view)", latest_samples_translated.sort('workdir_cluster_id'))

    # Check for situations like this, where A and B generated cluster IDs 000015 and 000033 respectively, but nothing from 000033
    # ended up in the persistent script.
    # This probably won't happen anymore due to better handling of decimated clusters, buuuuuuut just in case...
    #
    # samp | dis | workdir | last? | cluster_id
    # -----------------------------------------
    #  A   |  5  | 000015  | true  | 000033
    #  B   |  5  | 000033  | false | 000033
    #
    multi_workdir_newnames = (latest_samples_translated.group_by("cluster_id")
        .agg(pl.col("workdir_cluster_id").n_unique().alias("n_workdirs"))
        .filter(pl.col("n_workdirs") > 1)
        .get_column("cluster_id")
        .to_list()
    )

    # Must be in multi_workdir_newnames too, or else we will break persistent cluster IDs
    # that are actually working properly
    mask = (
        latest_samples_translated["cluster_id"].is_in(multi_workdir_newnames) &
        (latest_samples_translated["workdir_cluster_id"] == latest_samples_translated["cluster_id"])
    )

    problematic_stuff = latest_samples_translated.filter(mask)
    non_problematic_stuff = latest_samples_translated.filter(~mask) # we will add it back in later!!
    print_df_to_debug_log("Problematic stuff", problematic_stuff)

    # We want each problematic workdir_cluster_id to be given a new cluster_id (ie, a new persistent cluster ID). Currently we might have this,
    # where DDDDDDDD and EEEEEEEE were in a 20 cluster last run that ended up getting split in this run, resulting in them keeping their latest_id,
    # which is already in use elsewhere.
    #
    # ┌───────────┬──────────────────┬────────────────────┬────────────────────────┬────────────────────────┬───────────────────────┬────────────┐
    # │ sample_id ┆ cluster_distance ┆ workdir_cluster_id ┆ in_20_cluster_last_run ┆ in_10_cluster_last_run ┆ in_5_cluster_last_run ┆ cluster_id │
    # │ ---       ┆ ---              ┆ ---                ┆ ---                    ┆ ---                    ┆ ---                   ┆ ---        │
    # │ str       ┆ i64              ┆ str                ┆ bool                   ┆ bool                   ┆ bool                  ┆ str        │
    # ╞═══════════╪══════════════════╪════════════════════╪════════════════════════╪════════════════════════╪═══════════════════════╪════════════╡
    # │ AAAAAAAA  ┆ 5                ┆ 000033             ┆ null                   ┆ null                   ┆ false                 ┆ 000033     │
    # │ BBBBBBBB  ┆ 5                ┆ 000033             ┆ null                   ┆ null                   ┆ false                 ┆ 000033     │
    # │ CCCCCCCC  ┆ 5                ┆ 000033             ┆ null                   ┆ null                   ┆ false                 ┆ 000033     │
    # │ DDDDDDDD  ┆ 20               ┆ 000104             ┆ true                   ┆ null                   ┆ null                  ┆ 000104     │
    # │ EEEEEEEE  ┆ 20               ┆ 000104             ┆ true                   ┆ null                   ┆ null                  ┆ 000104     │
    # └───────────┴──────────────────┴────────────────────┴────────────────────────┴────────────────────────┴───────────────────────┴────────────┘
    #
    # This would just give everything the same cluster_id (where next_cluster_id is zfilled max_cluster_id_as_int() + 1)
    #latest_samples_translated = latest_samples_translated.with_columns([
    #    pl.when(mask)
    #    .then(pl.lit(next_cluster_id(latest_samples_translated)))
    #    .otherwise(pl.col("cluster_id"))
    #    .alias("cluster_id")
    #])
    #
    # So we'll instead use... another group!
    if problematic_stuff.shape[0] > 0:
        problematic_stuff_grouped = problematic_stuff.group_by("cluster_id").agg(
            pl.col("sample_id"),
            pl.col("workdir_cluster_id").unique(),
            pl.col("cluster_distance").unique(),
            pl.col("in_20_cluster_last_run"),
            pl.col("in_10_cluster_last_run"),
            pl.col("in_5_cluster_last_run")
        )
        # we do NOT want to use polars expressions to update cluster_id because we want it to be different for every grouped row
        # polars dataframes are immutable in this context though so we have to use another dataframe
        new_rows = []
        max_cluster_int = max_cluster_id_as_int(latest_samples_translated)
        for row in problematic_stuff_grouped.iter_rows(named=True):
            max_cluster_int += 1 
            row["cluster_id"] = str(max_cluster_int).zfill(6)
            new_rows.append(row)
        problematic_stuff_grouped = pl.DataFrame(new_rows)
        kaboom = problematic_stuff_grouped.explode(["sample_id", "in_20_cluster_last_run", "in_10_cluster_last_run", "in_5_cluster_last_run"])
        kaboom = kaboom.with_columns([
            kaboom["workdir_cluster_id"].list.get(0).alias("workdir_cluster_id"),
            kaboom["cluster_distance"].list.get(0).alias("cluster_distance")
        ])
        print_df_to_debug_log("problematic stuff after fixing", kaboom)
        common_cols = [col for col in non_problematic_stuff.columns if col in kaboom.columns]
        if "special_handling" in non_problematic_stuff.columns:
            kaboom = kaboom.with_columns(special_handling=pl.lit("silliness"))
            common_cols.append("special_handling")
            logging.debug("common_cols %s", common_cols)
        try:
            latest_samples_translated = pl.concat([non_problematic_stuff.select(common_cols), kaboom.select(common_cols)], how='align_full')
            # NOTE: we can get away with align_full here because non_problematic_stuff is mutually exclusive to problematic_stuff (from which comes kaboom)
            # any repeated sample IDs could cause problems!!
            print_df_to_debug_log("latest_samples_translated after accounting for weirdness (sorted by workdir_cluster_id in this view)", latest_samples_translated.sort('workdir_cluster_id'))
        except Exception: # pylint: disable=broad-exception-caught
            print("Encontered error trying to merge dataframes. Will print debug information then exit.")
            print(f"common_cols: {common_cols}")
            print_df_to_debug_log("kaboom.select(common_cols)", kaboom.select(common_cols))
            print_df_to_debug_log("non_problematic_stuff.select(common_cols)", non_problematic_stuff.select(common_cols))
            exit(999)

    # group by persistent cluster ID
    grouped = latest_samples_translated.group_by("cluster_id").agg(
        pl.col("sample_id"),
        pl.col("cluster_distance").n_unique().alias("distance_nunique"),
        pl.col("cluster_distance").unique().alias("distance_values"),
        pl.col("in_20_cluster_last_run").unique(),
        pl.col("in_10_cluster_last_run").unique(),
        pl.col("in_5_cluster_last_run").unique(),
        pl.col("workdir_cluster_id").n_unique().alias("worky-dirky_cluster_id_nunique"),
        pl.col("workdir_cluster_id").unique().alias("worky-dirky")
    )
    if (grouped["distance_nunique"] > 1).any():
        logging.error("Panic! Here is the grouped dataframe:")
        logging.error(grouped)
        raise ValueError("Some clusters have multiple unique cluster_distance values.")
    if (grouped["worky-dirky_cluster_id_nunique"] > 1).any():
        logging.error("Panic! Here is the grouped dataframe:") # panic because this should be handled earlier
        logging.error(grouped)
        raise ValueError("Some clusters have multiple unique workdir_cluster_id values.")

    grouped = grouped.with_columns(
        grouped["distance_values"].list.get(0).alias("cluster_distance")
    ).drop(["distance_nunique", "distance_values"])

    grouped = grouped.with_columns(
        grouped["worky-dirky"].list.get(0).alias("workdir_cluster_id")
    ).drop(["worky-dirky_cluster_id_nunique", "worky-dirky"]).sort('cluster_id')

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
    grouped = grouped.drop(['in_20_cluster_last_run', 'in_10_cluster_last_run', 'in_5_cluster_last_run']) # will be readded upon join
    print_df_to_debug_log("After grouping and then intager-a-fy", grouped)

    grouped = grouped.drop("sample_id") # prevent creation of sample_id_right, also this is redundant when we agg() again later
    hella_redundant = (latest_samples_translated.drop("cluster_distance")).join(grouped, on="cluster_id")
    assert_series_equal(hella_redundant.select("workdir_cluster_id").to_series(), hella_redundant.select("workdir_cluster_id_right").to_series(), check_names=False)
    hella_redundant = hella_redundant.drop("workdir_cluster_id_right")
    grouped = None
    hella_redundant = hella_redundant.with_columns(
        pl.lit(None).cast(pl.Utf8).alias("cluster_parent"),
        pl.lit([]).cast(pl.List(pl.Utf8)).alias("cluster_children")
    )
    hella_redundant = hella_redundant.sort(["cluster_distance", "cluster_id"])
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
    #logging.debug("Sample map:")
    #logging.debug(sample_map)
    updates = []
    for row in hella_redundant.iter_rows(named=True):
        cluster_id, one_sample, distance = row["cluster_id"], row["sample_id"], row["cluster_distance"]
        logging.info("[{%s}] sample %s in cluster_id %s", distance, one_sample, cluster_id)
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
    #print_df_to_debug_log("updates", updates)
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
    #print_df_to_debug_log("after linking parents and children", hella_redundant)

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

        # this particular row is PER SAMPLE, not PER ClUSTER
        ~pl.coalesce(["in_20_cluster_last_run", "in_10_cluster_last_run", "in_5_cluster_last_run"]).alias("sample_newly_clustered")
    ]).drop("samples_previously_in_cluster")
    hella_redundant = hella_redundant.drop(["in_20_cluster_last_run", "in_10_cluster_last_run", "in_5_cluster_last_run"]).sort("cluster_id")

    # now let's get information as to which samples are new or old so we can highlight them
    hella_redundant = hella_redundant.with_columns(
        pl.col("sample_id").is_in(all_persistent_samples["sample_id"]).not_().alias("brand_new_sample")
    )
    print_df_to_debug_log("after processing what clusters and samples are brand new, sorted by cluster_id", hella_redundant)
    sample_level_information = hella_redundant.select(["sample_id", "cluster_distance", "cluster_id", "cluster_brand_new", "sample_newly_clustered", "brand_new_sample"])
    sample_level_information.write_csv(f'new_samples{today.isoformat()}.tsv', separator='\t')
    sample_level_information = None

    # We're grouping again! Is there a way to do this in the previous group? Maybe, but I'm trying
    # to get this up and running ASAP so whatever works, works!
    second_group = hella_redundant.group_by("cluster_id").agg(
        pl.col("workdir_cluster_id").unique().first(),
        pl.col("sample_id"),
        pl.col("cluster_distance").unique(),
        pl.col("cluster_brand_new").unique(),
        pl.col("cluster_needs_updating").unique(),
        pl.col("cluster_parent").unique(),
        pl.col("cluster_children").flatten().unique(),
        pl.col("brand_new_sample").any().alias("has_new_samples") # should be redundant with cluster_brand_new
    )

    # check cluster distances
    # TODO: this and other asserts will probably need to change if we change how we handle unclustered samples
    assert (second_group["cluster_distance"].list.len() == 1).all(), "cluster_distance lists have length ≠ 1"
    second_group = second_group.with_columns(pl.col("cluster_distance").list.get(0).alias("cluster_distance_int"))
    second_group = second_group.drop("cluster_distance").rename({"cluster_distance_int": "cluster_distance"})

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
    # TODO: why the hell was this ⬇️ a thing?
    #assert_series_equal(second_group.select("cluster_brand_new").to_series(), second_group.select("has_new_samples").to_series(), check_names=False)
    second_group.select(["cluster_id", "cluster_distance", "has_new_samples"]).write_csv(f"clusters_with_new_samples{today.isoformat()}.tsv", separator='\t')
    second_group = second_group.drop("has_new_samples")

    # convert [null] to null
    second_group = second_group.with_columns([
        pl.when(pl.col("cluster_children").list.get(0).is_null())
        .then(None)
        .otherwise(pl.col("cluster_children"))
        .alias("cluster_children"),

        pl.when(pl.col("cluster_brand_new"))
        .then(pl.lit(today.isoformat()))
        .otherwise(None)
        .alias("first_found")
    ])

    print_df_to_debug_log("after grouping hella_redundant by cluster_id, converting [null] to null, and other checks", second_group)

    # join with persistent cluster *metadata* tsv, which does not have lists of samples, since we already defined what cluster should
    # have what samples with this run. But for easier comparison, we'll also join our old friend persis_groupby_cluster.
    # TODO: eventually latest cluster metadata file should be joined here too
    persistent_clusters_meta = persistent_clusters_meta.with_columns(pl.lit(False).alias("cluster_brand_new"))
    all_cluster_information = second_group.join(persistent_clusters_meta, how="full", on="cluster_id", coalesce=True)
    all_cluster_information = all_cluster_information.with_columns([
        pl.when(pl.col("cluster_brand_new").is_null())
        .then(pl.when(pl.col("cluster_brand_new_right").is_null())
            .then(pl.lit(False))
            .otherwise(pl.col("cluster_brand_new_right"))
        )
        .otherwise(pl.col("cluster_brand_new"))
        .alias("cluster_brand_new"),

        pl.when(pl.col("first_found").is_null())
        .then(pl.col("first_found_right"))
        .otherwise(pl.col("first_found"))
        .alias("first_found"),

        ]).drop(["cluster_brand_new_right", "first_found_right"])
    all_cluster_information = all_cluster_information.join(persis_groupby_cluster, how="full", on="cluster_id", coalesce=True)
    all_cluster_information = all_cluster_information.with_columns(
         pl.when(pl.col("cluster_distance").is_null())
        .then(pl.when(pl.col("cluster_distance_right").is_null())
            .then(pl.lit("ERROR!!!!")) # should never happen
            .otherwise(pl.col("cluster_distance_right"))
        )
        .otherwise(pl.col("cluster_distance"))
        .alias("cluster_distance")
    ).drop("cluster_distance_right")

    all_cluster_information = all_cluster_information.rename({"sample_id_right": "sample_id_previously"})
    all_cluster_information = get_nwk_and_matrix_plus_local_mask(all_cluster_information, args.combineddiff).sort("cluster_id")

    if "last_update" not in all_cluster_information.columns:
        all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("last_update"))
    if "first_found" not in all_cluster_information.columns:
        all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("first_found"))

    # okay, everything looks good so far. let's get some URLs!!
    # we already asserted that token is defined with yes_microreact hence possibly-used-before-assignment can be turned off there
    if args.yes_microreact:
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
                #all_cluster_information = update_first_found(all_cluster_information, this_cluster_id) # already did that earlier
                all_cluster_information = update_last_update(all_cluster_information, this_cluster_id)
                if distance == 20 and not has_children:
                    logging.info("%s is a brand-new 20-cluster with no children, not uploading", this_cluster_id)
                    #all_cluster_information = update_first_found(all_cluster_information, this_cluster_id) # already did that earlier
                    all_cluster_information = update_last_update(all_cluster_information, this_cluster_id)
                    all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "cluster_needs_updating", False)
                    continue
                URL = create_new_mr_project(token, this_cluster_id) # pylint: disable=possibly-used-before-assignment
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
                    logging.debug("%s is a %s-cluster with no new samples (ergo no new children), skipping", this_cluster_id, distance)
                else:
                    logging.warning("%s has no samples! This is likely a decimated cluster that lost all of its samples. We will not be updating its MR project.", this_cluster_id) # TODO but we should
                    logging.warning("Row information:")
                    logging.warning(row)
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
                first_found_shorthand = f'{str(first_found.year)}{str(first_found.strftime("%b")).zfill(2)}'
            except TypeError: # fires if !today and first_found is datetime.date (as opposed to str)
                first_found = today if row["first_found"] is None else str(row["first_found"].isoformat())
                first_found_shorthand = f'{str(row["first_found"].year)}{str(row["first_found"].strftime("%b")).zfill(2)}'

            # Because there is never a situation where a new child cluster pops up in a parent cluster that doesn't need to be updated,
            # and because MR URLs don't need to be updated, clusters that don't need updating don't need to know parent/child URLs.
            if not needs_updating:
                if row["sample_id"] is not None:
                    logging.info("%s is a %s-cluster marked as not needing updating (no new samples and/or childless 20-cluster), skipping", this_cluster_id, distance)
                else:
                    logging.warning("%s has no samples! This is likely a decimated cluster that lost all of its samples. We will not be updating its MR project.", this_cluster_id) # TODO but we should
                continue
            
            with open("./REALER_template.json", "r") as real_template_json:
                mr_document = json.load(real_template_json)

            # project title
            fullID = f"{first_found_shorthand}-{str(distance).zfill(2)}SNP-{this_cluster_id}"
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
            update_existing_mr_project(token, URL, mr_document, 0)
            if args.shareemail is not None:
                share_mr_project(token, URL, args.shareemail) 

        all_cluster_information = all_cluster_information.sort("cluster_id")
        logging.info("Finished. All cluster information:")
        logging.info(all_cluster_information)
        os.remove(args.persistentclustermeta)
        os.remove(args.persistentids)
        os.remove(args.token)
        all_cluster_information.write_ndjson(f'all_cluster_information{today.isoformat()}.json')
        new_persistent_meta = all_cluster_information.select(['cluster_id', 'first_found', 'last_update', 'jurisdictions', 'microreact_url'])
    else:
        all_cluster_information = all_cluster_information.sort("cluster_id")
        logging.info("Not touching Microreact. Here's the final data table:")
        logging.info(all_cluster_information)
        all_cluster_information.write_ndjson(f'all_cluster_information{today.isoformat()}.json')
        new_persistent_meta = all_cluster_information.select(['cluster_id', 'first_found', 'last_update', 'jurisdictions'])

    # regardless of MR, these get written
    new_persistent_meta.write_csv(f'persistentMETA{today.isoformat()}.tsv', separator='\t')
    new_persistent_ids = hella_redundant.select(['sample_id', 'cluster_id', 'cluster_distance'])
    new_persistent_ids.write_csv(f'persistentIDS{today.isoformat()}.tsv', separator='\t')
    samp_persistent20cluster = new_persistent_ids.filter(pl.col('cluster_distance') == 20).select(['sample_id', 'cluster_id'])
    samp_persistent20cluster.write_csv(f'samp_persiscluster{today.isoformat()}.tsv', separator='\t')
    change_report = []

    for row in all_cluster_information.iter_rows(named=True):
        logging.debug("Checking %s", row['cluster_id'])
        try:
            what_is = set(row["sample_id"])
        except TypeError:
            what_is = set()
        try:
            what_was = set(row["sample_id_previously"])
        except TypeError:
            what_was = set()
        gained = list(what_is - what_was)
        lost = list(what_was - what_is)
        logging.debug("what is: %s", what_is)
        logging.debug("what was: %s", what_was)
        logging.debug("gained: %s", gained)
        logging.debug("lost: %s", lost)
        change_report.append({"cluster": f"{row['cluster_id']}@{row['cluster_distance']}", "gained": gained, "lost": lost, "maintained": list(what_is.intersection(what_was))})
    change_report_df = pl.DataFrame(change_report).with_columns([
        pl.when(pl.col('gained').list.len() == 0).then(None).otherwise(pl.col('gained')).alias("gained"),
        pl.when(pl.col('lost').list.len() == 0).then(None).otherwise(pl.col('lost')).alias("lost"),
        pl.when(pl.col('maintained').list.len() == 0).then(None).otherwise(pl.col('maintained')).alias("maintained"),
        ])
    logging.info("Finished. Here's how clusters have changed:")
    print(change_report_df)
    change_report_df.write_ndjson(f'change_report{today.isoformat()}.json')


def add_col_if_not_there(dataframe, column):
    if column not in dataframe.columns:
        return dataframe.with_columns(pl.lit(None).alias(column))
    return dataframe

def get_nwk_and_matrix_plus_local_mask(big_ol_dataframe, combineddiff):
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "a_matrix")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "a_tree")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_matrix")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_tree")
    #print_df_to_debug_log("get_nwk_and_matrix_plus_local_mask() got this dataframe", big_ol_dataframe)
    for row in big_ol_dataframe.iter_rows(named=True):
        this_cluster_id = row["cluster_id"]
        workdir_cluster_id = row["workdir_cluster_id"]
        if workdir_cluster_id is not None:
            amatrix = f"a{workdir_cluster_id}_dmtrx.tsv" if os.path.exists(f"a{workdir_cluster_id}_dmtrx.tsv") else None
            atree = f"a{workdir_cluster_id}.nwk" if os.path.exists(f"a{workdir_cluster_id}.nwk") else None
            if workdir_cluster_id != this_cluster_id: # do NOT remove this check
                if amatrix is not None:
                    if not os.path.exists(f"a{this_cluster_id}_dmtrx.tsv"): # and that's why we can't remove aforementioned check
                        os.rename(f"a{workdir_cluster_id}_dmtrx.tsv", f"a{this_cluster_id}_dmtrx.tsv")
                        amatrix = f"a{this_cluster_id}_dmtrx.tsv"
                        logging.debug("[%s] a_matrix was a%s_dmtrx.tsv, now a%s_dmtrx.tsv", this_cluster_id, workdir_cluster_id, this_cluster_id)
                    else:
                        logging.warning("[%s] Cannot rename a%s_dmtrx.tsv to a%s_dmtrx.tsv as that already exists; will maintain workdir name", this_cluster_id, workdir_cluster_id, this_cluster_id)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "a_matrix", amatrix)
                else:
                    logging.warning("[%s] workdir_cluster_id is %s but could not find a%s.nwk", this_cluster_id, workdir_cluster_id, workdir_cluster_id)

                if atree is not None:
                    if not os.path.exists(f"a{this_cluster_id}.nwk"):
                        os.rename(f"a{workdir_cluster_id}.nwk", f"a{this_cluster_id}.nwk")
                        atree = f"a{this_cluster_id}.nwk"
                        logging.debug("[%s] a_tree was a%s.nwk, now a%s.nwk", this_cluster_id, workdir_cluster_id, this_cluster_id)
                    else:
                       logging.warning("[%s] Cannot rename a%s.nwk to a%s.nwk as that already exists; will maintain workdir name", this_cluster_id, workdir_cluster_id, this_cluster_id)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "a_tree", atree)
                else:
                    logging.warning("[%s] workdir_cluster_id is %s but could not find a%s.nwk", this_cluster_id, workdir_cluster_id, workdir_cluster_id)
            else:
                logging.debug("[%s] assigned a_matrix and a_tree (workdir id matches cluster id)", this_cluster_id)
                big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "a_matrix", amatrix)
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
                        subprocess.run(f"python3 /scripts/find_clusters.py {btreepb} --type BM --collection-name {this_cluster_id} -jmatsu", shell=True, check=True)
                        logging.debug("[%s] ran find_clusters.py, looks like it returned 0", this_cluster_id)
                        bmatrix = f"b{this_cluster_id}_dmtrx.tsv" if os.path.exists(f"b{this_cluster_id}_dmtrx.tsv") else None
                    except subprocess.CalledProcessError as e:
                        logging.warning("[%s] Failed to generate locally-masked tree/matrix: %s", this_cluster_id, e.output)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_matrix", bmatrix)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_tree", btree)
        else:
            logging.debug("[%s] No workdir_cluster_id, this is probably a decimated cluster", this_cluster_id)
    print_df_to_debug_log("returning this dataframe from get_nwk_and_matrix_plus_local_mask()", big_ol_dataframe)
    return big_ol_dataframe

def update_cluster_column(df, cluster_id, column, new_value):
    assert column in df.columns, f"Tried to update {column} with {new_value} but {column} not in dataframe?"
    return df.with_columns(
        pl.when(df["cluster_id"] == cluster_id)
        .then(pl.lit(new_value))
        .otherwise(df[column])
        .alias(column))

def update_first_found(df: pl.DataFrame, cluster_id: str) -> pl.DataFrame:
    return df.with_columns(
        pl.when(df["cluster_id"] == cluster_id)
        .then(pl.lit(today.isoformat()))
        .otherwise(df["first_found"])
        .alias("first_found"))

def update_last_update(df: pl.DataFrame, cluster_id: str) -> pl.DataFrame:
    return df.with_columns(
        pl.when(df["cluster_id"] == cluster_id)
        .then(pl.lit(today.isoformat()))
        .otherwise(df["last_update"])
        .alias("last_update"))

def print_df_to_debug_log(dataframe_name: str, actual_dataframe: pl.DataFrame) -> None:
    logging.debug("%s", dataframe_name)
    logging.debug(actual_dataframe)

def update_existing_mr_project(token, mr_url, mr_document, retries=-1):
    retries += 1
    if retries < 3:
        update_resp = requests.post("https://microreact.org/api/projects/update",
            headers={"Access-Token": token, "Content-Type": "application/json; charset=UTF-8"},
            params={"project": mr_url, "access": "private"},
            timeout=100,
            json=mr_document)
        if update_resp.status_code == 200:
            URL = update_resp.json()['id']
            logging.debug("Updated MR project with id %s", URL)
        else:
            logging.error("Failed to update MR project with id %s [code %s]: %s", mr_url, update_resp.status_code, update_resp.text)
            logging.error("Will retry...")
            update_existing_mr_project(token, mr_url, mr_document, retries)
    else:
        logging.error("Failed to update MR project with id %s after multiple retries. Something's broken.", mr_url)
        exit(1)

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

def get_atree_raw(cluster_name: str, big_ol_dataframe: pl.DataFrame):
    try:
        atree_series = big_ol_dataframe.filter(pl.col("cluster_id") == cluster_name).select("a_tree")
        atree = atree_series.item()
        with open(atree, "r") as nwk_file:
            return nwk_file.readline() # only need first line
    except (OSError, TypeError): # OSError: File Not Found, TypeError: None
        return "((INITIAL_SUBTREE_ERROR:1,REPORT_THIS_BUG_TO_ASH:1):1,DO_NOT_INCLUDE_PHI_IN_REPORT:1);"

def get_btree_raw(cluster_name: str, big_ol_dataframe: pl.DataFrame):
    try:
        btree_series = big_ol_dataframe.filter(pl.col("cluster_id") == cluster_name).select("b_tree")
        btree = btree_series.item()
        with open(btree, "r") as nwk_file:
            return nwk_file.readline() # only need first line
    except (OSError, TypeError):
        return "((MASKED_SUBTREE_ERROR:1,REPORT_THIS_BUG_TO_ASH:1):1,DO_NOT_INCLUDE_PHI_IN_REPORT:1);"

def nullfill_LR(polars_df: pl.DataFrame, left_col: str, right_col:str) -> pl.DataFrame:
    return polars_df.with_columns(pl.col(left_col).fill_null(pl.col(right_col))).drop(right_col)

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

def get_amatrix_raw(cluster_name: str, big_ol_dataframe: pl.DataFrame):
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

def get_bmatrix_raw(cluster_name: str, big_ol_dataframe: pl.DataFrame):
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

def get_cluster_ids_for_sample(df: pl.DataFrame, sample_id: str) -> list[str]:
    return (
        df.filter(pl.col("sample_id") == sample_id)
          .get_column("cluster_id")
          .unique()
          .to_list()
    )

def has_microreact_url(df: pl.DataFrame, cluster_id_val: str) -> bool:
    result = df.filter(pl.col("cluster_id") == cluster_id_val).select(
        (pl.col("microreact_url").is_not_null()).alias("has_url")
    )
    return result.item() if result.height > 0 else False

def get_other_samples_in_cluster(df: pl.DataFrame, cluster_id: str, exclude_sample_ids: list[str]) -> list[str]:
    return (
        df.filter(
            (pl.col("cluster_id") == cluster_id) &
            (~pl.col("sample_id").is_in(exclude_sample_ids))
        )
        .get_column("sample_id")
        .unique()
        .to_list()
    )

def max_cluster_id_as_int(df: pl.DataFrame) -> str:
    all_ids = df.select(
        pl.col("workdir_cluster_id").cast(pl.Utf8),
        pl.col("cluster_id").cast(pl.Utf8),
    ).unpivot().drop_nulls()

    max_val = (
        all_ids.select(
            pl.col("value")
            .filter(pl.col("value").str.len_chars() == 6)
            .filter(pl.col("value").str.contains(r"^\d+$"))
            .cast(pl.Int64)
            .max()
        )
        .item()
    )

    return max_val

if __name__ == "__main__":
    main()