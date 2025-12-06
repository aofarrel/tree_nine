VERSION = "0.3.12" # does not necessarily match Tree Nine git version
verbose = True   # set to False unless you can't dump the logs folder; be aware Terra's logger is very laggy
print(f"PROCESS CLUSTERS - VERSION {VERSION}")

# pylint: disable=too-many-statements,too-many-branches,simplifiable-if-expression,too-many-locals,too-complex,consider-using-tuple,broad-exception-caught
# pylint: disable=wrong-import-position,useless-suppression,multiple-statements,line-too-long,consider-using-sys-exit,duplicate-code

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
pl.Config.set_tbl_width_chars(200)
pl.Config.set_fmt_str_lengths(500)
pl.Config.set_fmt_table_cell_list_len(500)
today = datetime.utcnow().date() # I don't care if this runs past midnight, give everything the same day!
print(f"It's {today} in Thurles right now. Up Tipp!")
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("requests").setLevel(logging.WARNING)

if os.path.isfile("/scripts/marcs_incredible_script_update.pl"):
    script_path = "/scripts"
elif os.path.isfile("./scripts/marcs_incredible_script_update.pl"):
    script_path = "./scripts"
else:
    raise FileNotFoundError

def main():

    print("################# (1) INPUT HANDLING #################")
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
    parser.add_argument('--disable_decimated_failsafe', action='store_true', help='do not error if a cluster on MR becomes decimated')
    parser.add_argument('--no_cleanup', action='store_true', help="do not clean up input files (this may break delocalization on Terra; only use this for rapid debug runs)")

    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO)
    if args.persistentclustermeta and not args.persistentids:
        raise ValueError("You provided --persistentclustermeta but no --persistentids, you need both or neither")
    if args.persistentids and not args.persistentclustermeta:
        raise ValueError("You provided --persistentids but no --persistentclustermeta, you need both or neither")
    if not args.persistentids and not args.persistentclustermeta:
        start_over = True
        print("You have not provided persistent IDs nor persistent cluster metadata. This will restart clustering.")
    else:
        start_over = False

    all_latest_samples = pl.read_csv(args.latestsamples,
        separator="\t", 
        schema_overrides={"latest_cluster_id": pl.Utf8}).filter(pl.col("latest_cluster_id").is_not_null())

    if not start_over:
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
        with open(args.token, 'r', encoding="utf-8") as file:
            token = file.readline()
    debug_logging_handler_df("Loaded all_latest_samples", all_latest_samples, "input_all_latest_samples")
    if not start_over:
        debug_logging_handler_df("Loaded all_persistent_samples", all_persistent_samples, "input_all_persistent_samples")

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

    if not start_over:
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
        debug_logging_handler_txt("Set of all latest samples", "input_handling")
        debug_logging_handler_txt(all_latest_samples_set, "input_handling")
        debug_logging_handler_txt("Set of all persistent samples", "input_handling")
        debug_logging_handler_txt(all_persistent_samples_set, "input_handling")
        if all_persistent_samples_set.issubset(all_latest_samples_set):
            debug_logging_handler_txt("All persistent samples is a subset of all latest samples", "input_handling", 20)
        else:
            samples_missing_from_latest = all_persistent_samples_set - all_latest_samples_set # these are sets so this excludes samples exclusive to all_latest
            debug_logging_handler_txt(f"Samples appear to be missing from the latest run: {samples_missing_from_latest}", "input_handling", 30)
            if args.allsamples:
                all_input_samples_including_unclustered = args.allsamples.split(',')
            else:
                all_input_samples_including_unclustered = None
                debug_logging_handler_txt("Missing args.allsamples; can't be sure if missing samples are dropped because they no longer cluster or if they were never input.", "input_handling", 30)
            for sample in samples_missing_from_latest:
                if all_input_samples_including_unclustered is None:
                    pass
                elif sample in all_input_samples_including_unclustered:
                    debug_logging_handler_txt(f"{sample} is newly unclustered", "input_handling", 30)
                else:
                    debug_logging_handler_txt(f"{sample} seems to have been dropped from inputs", "input_handling", 30)
                # get persistent cluster ID regardless
                cluster_ids = get_cluster_ids_for_sample(all_persistent_samples, sample)
                for cluster in cluster_ids:
                    if len(get_other_samples_in_cluster(all_persistent_samples, cluster, samples_missing_from_latest)) <= 1:
                        # In theory we could handle this, in practice it's a massive pain in the neck and very easy to mess up!!
                        debug_logging_handler_txt(f"{cluster} is decimated thanks to losing all samples (or all but one)", "input_handling", 30)
                        # IF AND ONLY IF this is not on MR (which should only happen if this is a 20-cluster with no subclusters),
                        # we can live with this being a decimated cluster.
                        if not has_microreact_url(persistent_clusters_meta, cluster):
                            debug_logging_handler_txt(f"{cluster} already lacks a Microreact URL, so we can live with it being decimated", "input_handling", 20)
                        elif args.disable_decimated_failsafe:
                            debug_logging_handler_txt(f"{cluster} has an MR URL but we will accept it being decimated due to --disable_decimated_failsafe", "input_handling", 30)
                        else:
                            debug_logging_handler_txt(f"{cluster} has an MR URL and should never be decimated. Cannot continue.", "input_handling", 40)
                            exit(55)
                    else:
                        debug_logging_handler_txt(f"Dropped {sample} from {cluster} but that seems to be okay", "input_handling", 20)
        
        print("################# (2) 𓅀 𓁪 THE MARC PERRY ZONE 𓁫 𓀂 #################")
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
        debug_logging_handler_txt("Preparing to run the absolute legend's script...", "marc_perry", 20)
        filtered_latest_20 = all_latest_20.join(all_persistent_20.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
        filtered_latest_10 = all_latest_10.join(all_persistent_10.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
        filtered_latest_5 = all_latest_5.join(all_persistent_5.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
        filtered_persistent_20 = all_persistent_20.join(all_latest_20.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
        filtered_persistent_10 = all_persistent_10.join(all_latest_10.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
        filtered_persistent_5 = all_persistent_5.join(all_latest_5.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
        for distance, dataframe in {20: filtered_latest_20, 10: filtered_latest_10, 5: filtered_latest_5}.items():
            debug_logging_handler_df(f"Filtered latest {distance}", dataframe, "marc_perry")
        for distance, dataframe in {20: filtered_persistent_20, 10: filtered_persistent_10, 5: filtered_persistent_5}.items():
            debug_logging_handler_df(f"Filtered persistent {distance}", dataframe, "marc_perry")

        filtered_latest_20.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_20.tsv', separator='\t', include_header=False)
        filtered_latest_10.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_10.tsv', separator='\t', include_header=False)
        filtered_latest_5.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_5.tsv', separator='\t', include_header=False)
        filtered_persistent_20.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_20.tsv', separator='\t', include_header=False)
        filtered_persistent_10.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_10.tsv', separator='\t', include_header=False)
        filtered_persistent_5.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_5.tsv', separator='\t', include_header=False)

        debug_logging_handler_txt("Actually running scripts...", "marc_perry", 20)
        subprocess.run(f"perl {script_path}/marcs_incredible_script_update.pl filtered_persistent_20.tsv filtered_latest_20.tsv", shell=True, check=True, capture_output=True, text=True)
        subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_20.tsv", shell=True, check=True)
        subprocess.run(f"perl {script_path}/marcs_incredible_script_update.pl filtered_persistent_10.tsv filtered_latest_10.tsv", shell=True, check=True, capture_output=True, text=True)
        subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_10.tsv", shell=True, check=True)
        subprocess.run(f"perl {script_path}/marcs_incredible_script_update.pl filtered_persistent_5.tsv filtered_latest_5.tsv", shell=True, check=True, capture_output=True, text=True)
        subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_5.tsv", shell=True, check=True)

        # TODO: why are were we not running equalize tabs except when logging is debug?

        # debug print basic rosetta stones
        if logging.root.level == logging.DEBUG:
            for rock in ['rosetta_stone_20.tsv', 'rosetta_stone_10.tsv', 'rosetta_stone_5.tsv']:
                with open(rock, 'r', encoding="utf-8") as file:
                    debug_logging_handler_txt(f"---------------------\nContents of {rock} (before strip_tsv and equalize_tabs):\n", "marc_perry", 10)
                    debug_logging_handler_txt(list(file), "marc_perry", 10)
                    #subprocess.run(f"/bin/bash {script_path}/equalize_tabs.sh {rock}", shell=True, check=True)
                    
        # get more information about merges... if we have any!
        rock_pairs = {'rosetta_stone_20.tsv':'rosetta_stone_20_merges.tsv', 
                    'rosetta_stone_10.tsv':'rosetta_stone_10_merges.tsv', 
                    'rosetta_stone_5.tsv':'rosetta_stone_5_merges.tsv'}
        for rock, merge_rock in rock_pairs.items():
            if os.path.isfile(merge_rock):
                debug_logging_handler_txt(f"Found {merge_rock}, indicating clusters merged at this distance", "marc_perry", 20)
                debug_logging_handler_txt(f"---------------------\nContents of {merge_rock} (before strip_tsv and equalize_tabs):\n", "marc_perry", 10)
                with open(merge_rock, 'r', encoding="utf-8") as file:
                    debug_logging_handler_txt(list(merge_rock), "marc_perry", 10)
                    #subprocess.run(f"/bin/bash {script_path}/equalize_tabs.sh {rock}", shell=True, check=True)
                    subprocess.run(f"/bin/bash {script_path}/strip_tsv.sh {rock} {merge_rock}", shell=True, check=True)
            else:
                debug_logging_handler_txt(f"Did not find {merge_rock}, indicating clusters didn't merge at this distance", "marc_perry", 20)

        # we need schema_overrides or else cluster IDs can become non-zfilled i64
        # For some godforesaken reason, some versions of polars will throw `polars.exceptions.ComputeError: found more fields than defined in 'Schema'` even if we set
        # infer_schema = True with a hella large infer_schema_length. Idk why because the exact same file works perfectly fine on my local installation of polars (polars==1.27.0)
        # without even needing to set anything with infer_schema!! Not even a try-except with the except having a three column schema works!! Ugh!!!
        # TODO: is this because the docker is polars==1.26.0?
        # ---> WORKAROUND: equalize_tabs.sh
        debug_logging_handler_txt("Processing perl outputs...", "marc_perry", 20)
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

        debug_logging_handler_df("latest_samples_translated after pl.coalesce and check (sorted by workdir_cluster_id in this view)", 
            latest_samples_translated.sort('workdir_cluster_id'), "marc_perry")

        print("################# (3) SPECIAL HANDLING #################")

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
        debug_logging_handler_txt(f"multi_workdir_newnames list: {multi_workdir_newnames}", "special_handling", 10)

        # Must be in multi_workdir_newnames too, or else we will break persistent cluster IDs
        # that are actually working properly
        mask = (
            latest_samples_translated["cluster_id"].is_in(multi_workdir_newnames) &
            (latest_samples_translated["workdir_cluster_id"] == latest_samples_translated["cluster_id"])
        )

        problematic_stuff = latest_samples_translated.filter(mask)
        non_problematic_stuff = latest_samples_translated.filter(~mask) # we will add it back in later!!
        debug_logging_handler_df("Problematic stuff (cluster_id is in multi_workdir_newnames and workdir_cluster_id = cluster_id)", problematic_stuff, "special_handling")

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
            debug_logging_handler_df("problematic stuff after fixing", kaboom, "special_handling")
            common_cols = [col for col in non_problematic_stuff.columns if col in kaboom.columns]
            if "special_handling" in non_problematic_stuff.columns:
                kaboom = kaboom.with_columns(special_handling=pl.lit("silliness"))
                common_cols.append("special_handling")
            try:
                # pl.concat will fail if kaboom has any columns that are fully null (this can happen!), so
                # we're going to explictly cast those columns as booleans to prevent errors
                kaboom = kaboom.select(
                    pl.col('sample_id'),
                    pl.col('cluster_distance'),
                    pl.col('workdir_cluster_id'),
                    pl.col('cluster_id'),
                    pl.col('special_handling'),
                    pl.col("in_20_cluster_last_run").cast(pl.Boolean).alias("in_20_cluster_last_run"),
                    pl.col("in_10_cluster_last_run").cast(pl.Boolean).alias("in_10_cluster_last_run"),
                    pl.col("in_5_cluster_last_run").cast(pl.Boolean).alias("in_5_cluster_last_run"),
                )

                print(kaboom) # DEBUG DELETE LATER

                debug_logging_handler_txt("Casted some columns as booleans (just in case)", "special_handling", 20) # later make lv 10 
                latest_samples_translated = pl.concat([non_problematic_stuff.select(common_cols), kaboom.select(common_cols)], how='align_full')
                debug_logging_handler_txt("If you're reading this, we successfully concatenated non_problematic_stuff with kaboom.", "special_handling", 20)
                # NOTE: we can get away with align_full here because non_problematic_stuff is mutually exclusive to problematic_stuff (from which comes kaboom)
                # any repeated sample IDs could cause problems!!
                debug_logging_handler_df("latest_samples_translated after accounting for weirdness (sorted by workdir_cluster_id in this view)", 
                    latest_samples_translated.sort('workdir_cluster_id'), "special_handling")
            except Exception as e:
                logging.basicConfig(level=logging.DEBUG) # because Terra may not delocalize the files we need, but we don't always want debug logging b/c it slows down Terra
                debug_logging_handler_txt("Encontered error trying to merge dataframes. Will print debug information then exit.", "special_handling", 40)
                debug_logging_handler_txt(f"Error seems to have been: {e}", "special_handling", 40)
                debug_logging_handler_txt(f"common_cols: {common_cols}", "special_handling", 40)
                debug_logging_handler_df("kaboom.select(common_cols)", kaboom.select(common_cols), "special_handling")
                debug_logging_handler_df("non_problematic_stuff.select(common_cols)", non_problematic_stuff.select(common_cols), "special_handling")
                debug_logging_handler_txt("Attempting to zip logs...", "special_handling", 40)
                try:
                    subprocess.run("zip -r logs.zip ./logs", shell=True, check=True)
                except Exception as eeeeeee:
                    debug_logging_handler_txt(f"Caught {eeeeeee} attempting to zip logs. It's just not our day.", "special_handling", 40)
                exit(231)

        print("################# (4) FIRST GROUP (by persistent cluster ID) #################")
        debug_logging_handler_txt("Grouping by persistent cluster ID", "first_group", 20)
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
            debug_logging_handler_txt("Fatal error: At least one row has a value greater than 1 in column 'distance_nunique' ", "first_group", 40)
            debug_logging_handler_df("Grouped dataframe", grouped, "first_group")
            raise ValueError("Some clusters have multiple unique cluster_distance values.")
        if (grouped["worky-dirky_cluster_id_nunique"] > 1).any():
            # panic because this should be handled earlier
            debug_logging_handler_txt("Fatal error: At least one row has a value greater than 1 in column (deep sigh) 'worky-dirky_cluster_id_nunique'", "first_groupt", 40)
            debug_logging_handler_df("Grouped dataframe", grouped, "first_group")
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
        debug_logging_handler_df("After grouping and then intager-a-fy", grouped, "first_group")
        grouped = grouped.drop("sample_id")
        debug_logging_handler_txt("Dropped sample_id from grouped to prevent creation of sample_id_right (also it's redundant when we agg() again later)", "first_group", 10)
        hella_redundant = (latest_samples_translated.drop("cluster_distance")).join(grouped, on="cluster_id")
        debug_logging_handler_txt("Joined grouped with latest_samples_translated upon cluster_id to form hella_redundant", "first_group", 10)
        assert_series_equal(hella_redundant.select("workdir_cluster_id").to_series(), hella_redundant.select("workdir_cluster_id_right").to_series(), check_names=False)
        debug_logging_handler_txt("Asserted hella_redundant's workdir_cluster_id == hella_redundant's workdir_cluster_id_right", "first_group", 10)
        hella_redundant = hella_redundant.drop("workdir_cluster_id_right")
        grouped = None
        latest_samples_translated = None
        debug_logging_handler_txt("Dropped workdir_cluster_id_right from hella_redundant, cleared grouped variable, cleared latest_samples_translated variable", "first_group", 10)
    else:
        # Force all_latest_samples to look kind of like hella_redundant
        debug_logging_handler_txt("Skipped a ton of stuff, since we're cluster IDs starting over...", "first_group", 20)
        hella_redundant = all_latest_samples.with_columns([
            pl.col("latest_cluster_id").alias("workdir_cluster_id"),
            pl.col("latest_cluster_id").alias("cluster_id"),
            pl.lit("restart").alias("special_handling"),
            pl.lit([False]).alias("samples_previously_in_cluster"),
            pl.lit(False).alias("in_20_cluster_last_run"),
            pl.lit(False).alias("in_10_cluster_last_run"),
            pl.lit(False).alias("in_5_cluster_last_run")
        ])
    
    print("################# (5) LINK PARENTS AND CHILDREN #################")
    # this is a super goofy way to link parents and children, but it seems to work
    # basically, we're building a dataframe that has per-sample information, and a bunch of
    # duplicated per-cluster information. this is EXTREMELY redundant and in a normal world
    # we would not do this!
    hella_redundant = hella_redundant.with_columns(
        pl.lit(None).cast(pl.Utf8).alias("cluster_parent"),
        pl.lit([]).cast(pl.List(pl.Utf8)).alias("cluster_children")
    )
    hella_redundant = hella_redundant.sort(["cluster_distance", "cluster_id"])
    debug_logging_handler_txt("Prepared empty parent-child columns in hella_redundant and sorted by cluster_distance and cluster_id", "parental_guidance", 10)
    debug_logging_handler_df("hella_redundant at start of parental_guidance", hella_redundant, "parental_guidance")

    debug_logging_handler_txt("Linking parents and children...", "parental_guidance", 20)
    sample_map = {dist: {} for dist in [5, 10, 20]}
    for row in hella_redundant.iter_rows(named=True):
        sample_map[row["cluster_distance"]][row["sample_id"]] = row["cluster_id"]
    updates = []
    for row in hella_redundant.iter_rows(named=True):
        cluster_id, one_sample, distance = row["cluster_id"], row["sample_id"], row["cluster_distance"]
        debug_logging_handler_txt(f"[{distance}] {one_sample} in cluster {cluster_id}", "parental_guidance", 10)
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
    debug_logging_handler_txt("Generated updates list, now using it to update the dataframe...", "parental_guidance", 20)
    for cluster_id, col, value in updates:
        #debug_logging_handler_txt(f"For cluster {cluster_id}, col {col}, val {value} in updates", "parental_guidance", 10) # too verbose even for debug logging
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
    debug_logging_handler_df("hella_redundant after linking parents and children", hella_redundant, "parental_guidance")

    print("################# (6) RECOGNIZE (have I seen you before?) #################")
    debug_logging_handler_txt("Determining which clusters are brand new and/or need updating", "recognize", 20)
    hella_redundant = hella_redundant.with_columns([
        # When samples_previously_in_cluster is:
        # [True, False]/[False, True] --> some samples were in cluster previously     --> old cluster, needs updating
        # [False]                     --> no samples were in this cluster previously  --> new cluster, needs updating
        # [True]                      --> all samples were in this cluster previously --> old cluster, unchanged
        #
        # Strictly speaking you don't need to update something that didn't exist anymore but I gave this a good long think
        # and have decided that "needs updating" (ie pinging microreact) is functionally what matters here.
        #
        # TODO: "After grouping and then intager-a-fy" can report [false] which becomes [true] in hella_redundant if special_handling is silliness
        #
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
    if start_over:
        hella_redundant = hella_redundant.with_columns(pl.lit(True).alias("brand_new_sample"))
    else:
        hella_redundant = hella_redundant.with_columns(
            pl.col("sample_id").is_in(all_persistent_samples["sample_id"]).not_().alias("brand_new_sample")
        )
    debug_logging_handler_df("after processing what clusters and samples are brand new, sorted by cluster_id", hella_redundant, "recognize")
    sample_level_information = hella_redundant.select(["sample_id", "cluster_distance", "cluster_id", "cluster_brand_new", "sample_newly_clustered", "brand_new_sample"])
    sample_level_information.write_csv(f'all_samples{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote all_samples{today.isoformat()}.tsv from some of hella_redundant's columns", "recognize", 20)
    sample_level_information.filter(pl.col('brand_new_sample')).write_csv(f'new_samples{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote new_samples{today.isoformat()}.tsv which should only have the brand new samples in it", "recognize", 20)
    sample_level_information = None

    print("################# (7) SECOND GROUP (back at it again) #################")
    debug_logging_handler_txt("Grouping again...", "second_group", 20)
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
    debug_logging_handler_txt("Reformatting and performing checks...", "second_group", 20)
    assert (second_group["cluster_distance"].list.len() == 1).all(), "cluster_distance lists have length ≠ 1"
    debug_logging_handler_txt("Asserted all cluster_distance lists have a len of precisely 1", "second_group", 10)
    second_group = second_group.with_columns(pl.col("cluster_distance").list.get(0).alias("cluster_distance_int"))
    second_group = second_group.drop("cluster_distance").rename({"cluster_distance_int": "cluster_distance"})
    debug_logging_handler_txt("Converted lists of cluter distance (which we know have a len of 1) into ints", "second_group", 10)

    # check parenthood
    # in polars, [null] is considered to have a length of 1. I'm checking by distance rather than all clusters at once in case that changes.
    # be aware: https://github.com/pola-rs/polars/issues/18522
    debug_logging_handler_txt("Some checks rely upon len([pl.Null]) == 1, which is polars behavior that may change later. Beware!", "second_group", 30)
    assert ((second_group.filter(second_group["cluster_distance"] == 5))["cluster_parent"].list.len() == 1).all(), "5-cluster with multiple cluster_parents"
    assert ((second_group.filter(second_group["cluster_distance"] == 10))["cluster_parent"].list.len() == 1).all(), "10-cluster with multiple cluster_parents"
    assert ((second_group.filter(second_group["cluster_distance"] == 20))["cluster_parent"].list.len() == 1).all(), "20-cluster with cluster_parent *OR* preliminary null error"
    debug_logging_handler_txt("Asserted all clusters, if they have a parent, have only one parent", "second_group", 10)
    second_group = second_group.with_columns(pl.col("cluster_parent").list.get(0).alias("cluster_parent_str"))
    second_group = second_group.drop("cluster_parent").rename({"cluster_parent_str": "cluster_parent"})

    # check otherstuff
    assert ((second_group.filter(second_group["cluster_distance"] == 20))["cluster_parent"].is_null()).all(), "20-cluster with cluster_parent *OR* secondary null error"
    debug_logging_handler_txt("Asserted no 20 clusters have parents and that pl.Null hasn't exploded", "second_group", 10)
    assert ((second_group.filter(second_group["cluster_distance"] == 5))["cluster_children"].list.get(0).is_null()).all(), "5-cluster with cluster_children"
    debug_logging_handler_txt("Asserted no 5 clusters have children", "second_group", 10)
    assert (second_group["cluster_needs_updating"].list.len() == 1).all(), "Cluster not sure if it needs updating"
    debug_logging_handler_txt("Asserted all len(cluster_needs_updating) == 1, but this may not catch all edge cases involving cluster updating", "second_group", 10)
    second_group = second_group.with_columns(pl.col("cluster_needs_updating").list.get(0).alias("cluster_needs_updating_bool"))
    second_group = second_group.drop("cluster_needs_updating").rename({"cluster_needs_updating_bool": "cluster_needs_updating"})
    assert (second_group["cluster_brand_new"].list.len() == 1).all(), "Cluster not sure if it's new"
    debug_logging_handler_txt("Asserted all len(cluster_brand_new) == 1, but this may not catch all edge cases involving cluster newness", "second_group", 10)
    second_group = second_group.with_columns(pl.col("cluster_brand_new").list.get(0).alias("cluster_brand_new_bool"))
    second_group = second_group.drop("cluster_brand_new").rename({"cluster_brand_new_bool": "cluster_brand_new"})
    # TODO: why the hell was this ⬇️ a thing? was this to check it was NOT equal? ...should we maybe readd that actually?
    #assert_series_equal(second_group.select("cluster_brand_new").to_series(), second_group.select("has_new_samples").to_series(), check_names=False)
    second_group.select(["cluster_id", "cluster_distance", "has_new_samples"]).write_csv(f"clusters_with_new_samples{today.isoformat()}.tsv", separator='\t')
    second_group = second_group.drop("has_new_samples")

    # convert [null] to null
    debug_logging_handler_txt("Converting [null] to null...", "second_group", 20)
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

    debug_logging_handler_df("after grouping hella_redundant by cluster_id, converting [null] to null, and other checks", second_group, "second_group")

    print("################# (7) JOIN with persistent information #################")
    # join with persistent cluster *metadata* tsv, which does not have lists of samples, since we already defined what cluster should
    # have what samples with this run. But for easier comparison, we'll also join our old friend persis_groupby_cluster.
    # TODO: eventually latest cluster metadata file should be joined here too
    if start_over:
        debug_logging_handler_txt("Generating metadata fresh (since we're starting over)...", "join_metadata", 20)
        all_cluster_information = second_group.with_columns(cluster_brand_new=True, first_found=today)
    else:
        debug_logging_handler_txt("Joining with the persistent metadata TSV...", "join_metadata", 20)
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
        debug_logging_handler_df("after joining with persistent_clusters_meta", all_cluster_information, "join_metadata")

        # persis_groupby_cluster is only created if start_over is false, so it existing here is okay
        debug_logging_handler_txt("Joining persis_groupby_cluster...", "join_metadata", 20)
        all_cluster_information = all_cluster_information.join(persis_groupby_cluster, how="full", on="cluster_id", coalesce=True) # pylint: disable=possibly-used-before-assignment
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
        debug_logging_handler_df("after joining with persis_groupby_cluster and getting nwk, matrix, and mask", all_cluster_information, "join_metadata")

    if "last_update" not in all_cluster_information.columns:
        all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("last_update"))
    if "first_found" not in all_cluster_information.columns:
        all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("first_found"))
    if "jurisdictions" not in all_cluster_information.columns:
        all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("jurisdictions"))
    if "sample_id_previously" not in all_cluster_information.columns: # happens if start_over
        all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("sample_id_previously"))


    # okay, everything looks good so far. let's get some URLs!!
    # we already asserted that token is defined with yes_microreact hence possibly-used-before-assignment can be turned off there
    print("################# (8) MICROREACT #################")
    if args.yes_microreact:
        debug_logging_handler_txt("Assigning self-URLs...", "microreact", 20)
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
                    debug_logging_handler_txt(f"{this_cluster_id} is a brand-new 20-cluster with no children, not uploading", "microreact", 10)
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
                        debug_logging_handler_txt(f"{this_cluster_id} is an old 20-cluster with no children, not uploading", "microreact", 10)
                        all_cluster_information = update_last_update(all_cluster_information, this_cluster_id)
                        all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "cluster_needs_updating", False)
                        continue
                    debug_logging_handler_txt(f"{this_cluster_id} isn't brand new, but is flagged as needing an update and has no URL. Will make a new URL.", "microreact", 30)
                    URL = create_new_mr_project(token, this_cluster_id)
                    all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "microreact_url", URL)
                    all_cluster_information = update_last_update(all_cluster_information, this_cluster_id)
                else:
                    debug_logging_handler_txt(f"{this_cluster_id}'s URL ({URL}) seems valid, but we're not gonna check it", "microreact", 10)

            else:
                debug_logging_handler_txt(f"{this_cluster_id} seems unchanged", "microreact", 10)
                continue

        all_cluster_information = all_cluster_information.with_columns(
            pl.lit(None).alias("parent_URL"),
            pl.lit(None).alias("children_URLs") # TODO: how are gonna keep these ordered with the children column... does that even matter?
        )

        # now that everything has a URL, or doesn't need one, iterate a second time to get URLs of parents and children
        debug_logging_handler_txt("Searching for MR URLs of parents and children...", "microreact", 20)
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
                    debug_logging_handler_txt(f"{this_cluster_id}@{distance} has no new samples (ergo no new children), skipping", "microreact", 10)
                else:
                    debug_logging_handler_txt(f"{this_cluster_id}@{distance} has no samples! This is likely a decimated cluster that lost all of its samples. We will not be updating its MR project.", "microreact", 30)
                    debug_logging_handler_txt("Row information of this decimated cluster:", "microreact", 10)
                    debug_logging_handler_txt(row, "microreact", 10)
                    continue

            if has_parent:
                # get value of column "microreact_url" when column "cluster_id" matches the string in row["cluster_parent"]
                parent_URL = all_cluster_information.filter(
                    all_cluster_information["cluster_id"] == row["cluster_parent"]
                )["microreact_url"].item()
                all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "parent_URL", parent_URL)
                
            if has_children:
                # get value of column "microreact_url" when column "cluster_id" matches a string in the list in row["cluster_children"]
                children_URLs = all_cluster_information.filter(
                    all_cluster_information["cluster_id"].is_in(row["cluster_children"])
                )["microreact_url"].to_list()
                all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "children_URLs", children_URLs)
        parent_URL, children_URLs, URL = None, None, None
        debug_logging_handler_df("all_cluster_information prior to upload", all_cluster_information, "microreact")

        # now that everything can be crosslinked, iterate one more time to actually upload
        # yes three iterations is weird, cringe even. I'm doing this to make debugging easier.
        debug_logging_handler_txt("Splitting and uploading...", "microreact", 10)
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
                    debug_logging_handler_txt(f"{this_cluster_id}@{distance} marked as not needing updating (no new samples and/or childless 20-cluster), skipping", "microreact", 10)
                else:
                    debug_logging_handler_txt(f"You probably already know this, but {this_cluster_id}@{distance} has no samples!", "microreact", 30)
                continue
            
            with open("./REALER_template.json", "r", encoding="utf-8") as real_template_json:
                mr_document = json.load(real_template_json)

            # microreact needs all panals specified in a JSON file to be filled out, or else the workspace won't load

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

            # tree labels (metadata)
            if os.path.isfile("./metadata_combined.tsv"):
                # TODO: if MR cannot handle sample IDs being in the table that aren't on the tree, we will need to do more processing here
                debug_logging_handler_txt("Found metadata_combined.tsv, will use that for metadata", "microreact", 20)
                metadata_dict = csv.reader("./metadata_combined.tsv", delimiter="\t")
            else:
                debug_logging_handler_txt("Could not find metadata_combined.tsv, will mark as undefined per current CDPH guidelines", "microreact", 20)
                metadata_dict = [
                    {
                        "id": sample_id,
                        "Country": "UNDEFINED",
                        "Epi_Duplication": "UNDEFINED",
                        "Latitude": "UNDEFINED",
                        "Longitude": "UNDEFINED",
                        "Patient_County": "UNDEFINED",
                        "State": "UNDEFINED",
                        "Submitter_Facility": "UNDEFINED",
                        "Submitter_Facility_Sample_ID": "UNDEFINED",
                        "Year_Collected": "UNDEFINED",
                        "Sequencing_Facility": "UNDEFINED"
                    }
                    for sample_id in sample_id_list
                ]
            debug_logging_handler_txt(f"Metadata dictionary: {metadata_dict}", "microreact", 10)
            output = io.StringIO(newline='') # get rid of carriage return (this is kind of a silly way to do it but it works)
            writer = csv.DictWriter(output, fieldnames=metadata_dict[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(metadata_dict)
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

            debug_logging_handler_txt(f"MR document for {this_cluster_id}:", "microreact", 10)
            debug_logging_handler_txt(f"{mr_document}", "microreact", 30)

            # actually upload
            assert URL is not None, f"No Microreact URL for {this_cluster_id}!"
            update_existing_mr_project(token, URL, mr_document, 0)
            if args.shareemail is not None:
                share_mr_project(token, URL, args.shareemail) 

        all_cluster_information = all_cluster_information.sort("cluster_id")
        debug_logging_handler_txt("Finished uploading to Microreact", "microreact", 20)
        debug_logging_handler_df("all_cluster_information after MR uploads", all_cluster_information, "microreact")
        new_persistent_meta = all_cluster_information.select(['cluster_id', 'first_found', 'last_update', 'jurisdictions', 'microreact_url'])
    else:
        all_cluster_information = all_cluster_information.sort("cluster_id")
        debug_logging_handler_df("Not touching Microreact", all_cluster_information, "microreact")
        if "microreact_url" not in all_cluster_information.columns:
            all_cluster_information = all_cluster_information.with_columns(pl.lit(None).alias("microreact_url"))
        if "parent_URL" not in all_cluster_information.columns and "children_URLs" not in all_cluster_information.columns:
            all_cluster_information = all_cluster_information.with_columns(
                pl.lit(None).alias("parent_URL"),
                pl.lit(None).alias("children_URLs") # TODO: how are gonna keep these ordered with the children column... does that even matter?
            )
        new_persistent_meta = all_cluster_information.select(['cluster_id', 'first_found', 'last_update', 'jurisdictions'])


    print("################# (9) FINISHING UP #################")
    if not args.no_cleanup:
        if args.persistentclustermeta:
            os.remove(args.persistentclustermeta)
        if args.persistentids:
            os.remove(args.persistentids)
        if args.token is not None:
            os.remove(args.token)
        debug_logging_handler_txt("Deleted input persistentclustermeta, input persistentids, and input token", "final", 10)
    all_cluster_information.write_ndjson(f'all_cluster_information{today.isoformat()}.json')
    debug_logging_handler_txt(f"Wrote all_cluster_information{today.isoformat()}.json", "final", 20)
    
    # persistentMETA and persistentIDS are needed for subsequent runs
    new_persistent_meta.write_csv(f'persistentMETA{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote persistentMETA{today.isoformat()}.tsv", "final", 20)
    new_persistent_ids = hella_redundant.select(['sample_id', 'cluster_id', 'cluster_distance'])
    new_persistent_ids.write_csv(f'persistentIDS{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote persistentIDS{today.isoformat()}.tsv", "final", 20)

    # the sample \t cluster TSVs can be used to convert to annotated nextstrain format
    samp_persistent20cluster = new_persistent_ids.filter(pl.col('cluster_distance') == 20).select(['sample_id', 'cluster_id'])
    samp_persistent10cluster = new_persistent_ids.filter(pl.col('cluster_distance') == 10).select(['sample_id', 'cluster_id'])
    samp_persistent5cluster = new_persistent_ids.filter(pl.col('cluster_distance') == 5).select(['sample_id', 'cluster_id'])
    samp_persistent20cluster.write_csv(f'samp_persis20cluster{today.isoformat()}.tsv', separator='\t')
    samp_persistent10cluster.write_csv(f'samp_persis10cluster{today.isoformat()}.tsv', separator='\t')
    samp_persistent5cluster.write_csv(f'samp_persis5cluster{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote samp_persis20cluster{today.isoformat()}.tsv, samp_persis10cluster{today.isoformat()}.tsv, and samp_persis5cluster{today.isoformat()}.tsv", "final", 20)
    
    change_report = []
    debug_logging_handler_txt("Building change report...", "final", 20)
    for row in all_cluster_information.iter_rows(named=True):
        try:
            what_is = set(row["sample_id"])
        except TypeError:
            what_is = set()
        try:
            what_was = set(row["sample_id_previously"])
        except (TypeError, KeyError):
            what_was = set()
        gained = list(what_is - what_was)
        lost = list(what_was - what_is)
        change_report.append({"cluster": f"{row['cluster_id']}@{row['cluster_distance']}", 
            "gained": gained, "lost": lost, "kept": list(what_is.intersection(what_was)),
            "microreact_url": row['microreact_url'], "dist": row['cluster_distance']})
    change_report_df = pl.DataFrame(change_report).with_columns([
        pl.when(pl.col('gained').list.len() == 0).then(None).otherwise(pl.col('gained')).alias("gained"),
        pl.when(pl.col('lost').list.len() == 0).then(None).otherwise(pl.col('lost')).alias("lost"),
        pl.when(pl.col('kept').list.len() == 0).then(None).otherwise(pl.col('kept')).alias("kept"),
        pl.when(pl.col('microreact_url').is_null()).then(None).otherwise(pl.col('microreact_url')).alias("microreact_url"),
    ])
    change_report_df = change_report_df.with_columns([
        pl.when(pl.col('gained').is_not_null()).then(pl.col('gained').list.len()).otherwise(pl.lit(0)).alias("n_gained"),
        pl.when(pl.col('lost').is_not_null()).then(pl.col('lost').list.len()).otherwise(pl.lit(0)).alias("n_lost"),
        pl.when(pl.col('kept').is_not_null()).then(pl.col('kept').list.len()).otherwise(pl.lit(0)).alias("n_kept"),
    ])
    change_report_df = change_report_df.with_columns(n_now=pl.col('n_gained')-pl.col('n_lost')+pl.col('n_kept'))

    # sometimes it thinks the dist column is a string, let's cast it as an int
    change_report_df = change_report_df.with_columns(pl.col("dist").cast(pl.Int32))
    change_report_df_no_twenties = change_report_df.filter(pl.col('dist') != 20)

    pl.Config.set_tbl_width_chars(200)
    with open(f"change_report_full{today.isoformat()}.txt", "a", encoding="utf-8") as full:
        with open(f"change_report_cdph{today.isoformat()}.txt", "a", encoding="utf-8") as cdph:
            full.write("Existing clusters that lost samples (note: it's possible to gain and lose)\n")
            print(change_report_df.filter(pl.col("lost").is_not_null()).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'microreact_url', 'lost']), file=full)
            cdph.write("Existing clusters that lost samples (note: it's possible to gain and lose)\n")
            print(change_report_df_no_twenties.filter(pl.col("lost").is_not_null()).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'microreact_url', 'lost']), file=cdph)

            full.write("Existing clusters that gained samples (note: it's possible to gain and lose)\n")
            print(change_report_df.filter((pl.col("gained").is_not_null().and_(pl.col("kept").is_not_null()))).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'n_now', 'microreact_url', 'gained']), file=full)
            cdph.write("Existing clusters that gained samples (note: it's possible to gain and lose)\n")
            print(change_report_df_no_twenties.filter((pl.col("gained").is_not_null().and_(pl.col("kept").is_not_null()))).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'n_now', 'microreact_url', 'gained']), file=cdph)

            full.write("Brand new clusters\n")
            print(change_report_df.filter((pl.col("gained").is_not_null().and_(pl.col("kept").is_null()))).select(['cluster', 'n_gained', 'microreact_url', 'gained']), file=full)
            cdph.write("Brand new clusters\n")
            print(change_report_df_no_twenties.filter((pl.col("gained").is_not_null().and_(pl.col("kept").is_null()))).select(['cluster', 'n_gained', 'microreact_url', 'gained']), file=cdph)

            full.write("Decimated clusters\n")
            print(change_report_df.filter((pl.col("lost").is_null()).and_(pl.col("gained").is_null()).and_(pl.col("kept").is_null())).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'n_now', 'microreact_url']), file=full)
            cdph.write("Decimated clusters\n")
            print(change_report_df_no_twenties.filter((pl.col("lost").is_null()).and_(pl.col("gained").is_null()).and_(pl.col("kept").is_null())).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'n_now', 'microreact_url']), file=cdph)

            full.write("Unchanged clusters\n")
            print(change_report_df.filter((pl.col("lost").is_null()).and_(pl.col("gained").is_null()).and_(pl.col("kept").is_not_null())).select(['cluster', 'n_now', 'microreact_url']), file=full)
            cdph.write("Unchanged clusters\n")
            print(change_report_df.filter((pl.col("lost").is_null()).and_(pl.col("gained").is_null()).and_(pl.col("kept").is_not_null())).select(['cluster', 'n_now', 'microreact_url']), file=cdph)


    change_report_df.write_ndjson(f'change_report{today.isoformat()}.json')
    debug_logging_handler_txt(f"Finished. Saved change report dataframe as change_report{today.isoformat()}.json", "final", 20)

# Ultimately we want to have our cake and eat it too with regard to logging.
# 1) Terra might not delocalize log files if a task exits early
# 2) We want multiple log files when possible to make debugging actually feasible

def debug_logging_handler_txt(msg: str, logfile: str, loglevel=10):
    time = datetime.utcnow().strftime("%H:%M")
    if loglevel == 40:
        logging.error("[%s @ %s] %s", logfile, time, msg)
    elif loglevel == 30:
        logging.warning("[%s @ %s] %s", logfile, time, msg)
    elif loglevel == 20:
        logging.info("[%s @ %s] %s", logfile, time, msg)
    else:
        logging.debug("[%s @ %s] %s", logfile, time, msg)
    try:
        with open("./logs/"+logfile+".log", "a", encoding="utf-8") as f:
            f.write(str(msg) + "\n")
    except Exception:
        logging.warning("Logging error!")
        logging.warning(f"   msg = {msg}") # pylint: disable=logging-fstring-interpolation
        logging.warning(f"   logfile = {logfile}") # pylint: disable=logging-fstring-interpolation
        logging.warning(f"   loglevel = {loglevel}") # pylint: disable=logging-fstring-interpolation
        exit(1)

def debug_logging_handler_df(title: str, dataframe: pl.DataFrame, logfile: str):
    time = datetime.utcnow().strftime("%H:%M")
    logging.info("[%s @ %s] Dumping debug dataframe with info: %s", logfile, time, title)
    try:
        json_name = logfile+"__"+title.replace(" ", "_")
        dataframe.write_ndjson("./logs/"+json_name+'.json')
        logging.info("[%s @ %s] SEE ALSO: %s.json", logfile, time, json_name)
    except Exception: # ignore: broad-exception-caught
        logging.info("[%s @ %s] Failed to write json version of dataframe, rely on polars' best efforts below (this is a logging error and is probably fine)", logfile, time)
        with open("./logs/"+logfile+".log", "a", encoding="utf-8") as f:
            f.write(title + "\n")
            f.write(dataframe)
    # in case of early exit, ALSO dump to stderr if logging.debug
    logging.debug(dataframe)
    

def add_col_if_not_there(dataframe: pl.DataFrame, column: str):
    if column not in dataframe.columns:
        return dataframe.with_columns(pl.lit(None).alias(column))
    return dataframe

def get_nwk_and_matrix_plus_local_mask(big_ol_dataframe: pl.DataFrame, combineddiff: str):
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "a_matrix")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "a_tree")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_matrix")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_tree")
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
                        subprocess.run(f"python3 {script_path}/find_clusters.py {btreepb} --type BM --collection-name {this_cluster_id} -jmatsu", shell=True, check=True)
                        logging.debug("[%s] ran find_clusters.py, looks like it returned 0", this_cluster_id)
                        bmatrix = f"b{this_cluster_id}_dmtrx.tsv" if os.path.exists(f"b{this_cluster_id}_dmtrx.tsv") else None
                    except subprocess.CalledProcessError as e:
                        logging.warning("[%s] Failed to generate locally-masked tree/matrix: %s", this_cluster_id, e.output)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_matrix", bmatrix)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_tree", btree)
        else:
            logging.debug("[%s] No workdir_cluster_id, this is probably a decimated cluster", this_cluster_id)
    return big_ol_dataframe

def update_cluster_column(df: pl.DataFrame, cluster_id, column, new_value):
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
            debug_logging_handler_txt(f"Updated MR project {URL}", "microreact", 10)
        else:
            debug_logging_handler_txt(f"Failed to update MR project {mr_url} [code {update_resp.status_code}]: {update_resp.text}", "microreact", 30)
            debug_logging_handler_txt("RETRYING...", "microreact", 20)
            update_existing_mr_project(token, mr_url, mr_document, retries)
    else:
        debug_logging_handler_txt(f"Failed to update MR project {mr_url} after multiple retries. Something's broken.", "microreact", 40)
        exit(1)

def share_mr_project(token, mr_url, email):
    api_url = "https://microreact.org/api/shares/add"
    params = {"id": mr_url}
    headers = {"Access-Token": token}
    data = {"emails": [email], "role": "viewer" }
    share_resp = requests.post(api_url, headers=headers, params=params, json=data, timeout=100)
    if share_resp.status_code != 200: 
        debug_logging_handler_txt(f"Failed to share MR project {mr_url} [code {share_resp.status_code}]: {share_resp.text}", "microreact", 40)
        debug_logging_handler_txt("NOT retrying as this is probably a permissions issue.", "microreact", 40)

def create_new_mr_project(token, this_cluster_id):
    with open("./BLANK_template.json", "r", encoding="utf-8") as temp_proj_json:
        mr_document = json.load(temp_proj_json)
    update_resp = requests.post("https://microreact.org/api/projects/create",
        headers={"Access-Token": token, "Content-Type": "application/json; charset=UTF-8"},
        params={"access": "private"},
        timeout=100,
        json=mr_document)
    if update_resp.status_code == 200:
        URL = update_resp.json()['id']
        debug_logging_handler_txt(f"{this_cluster_id} is brand new and has been assigned {URL} and a first-found date", "microreact", 10)
        return URL
    debug_logging_handler_txt(f"Failed to create new MR project for {this_cluster_id} [code {update_resp.status_code}]: {update_resp.text}", "microreact", 40)
    return None

def get_atree_raw(cluster_name: str, big_ol_dataframe: pl.DataFrame):
    try:
        atree_series = big_ol_dataframe.filter(pl.col("cluster_id") == cluster_name).select("a_tree")
        atree = atree_series.item()
        with open(atree, "r", encoding="utf-8") as nwk_file:
            return nwk_file.readline() # only need first line
    except (OSError, TypeError): # OSError: File Not Found, TypeError: None
        return "((INITIAL_SUBTREE_ERROR:1,REPORT_THIS_BUG_TO_ASH:1):1,DO_NOT_INCLUDE_PHI_IN_REPORT:1);"

def get_btree_raw(cluster_name: str, big_ol_dataframe: pl.DataFrame):
    try:
        btree_series = big_ol_dataframe.filter(pl.col("cluster_id") == cluster_name).select("b_tree")
        btree = btree_series.item()
        with open(btree, "r", encoding="utf-8") as nwk_file:
            return nwk_file.readline() # only need first line
    except (OSError, TypeError):
        return "((MASKED_SUBTREE_ERROR:1,REPORT_THIS_BUG_TO_ASH:1):1,DO_NOT_INCLUDE_PHI_IN_REPORT:1);"

def nullfill_LR(polars_df: pl.DataFrame, left_col: str, right_col:str) -> pl.DataFrame:
    return polars_df.with_columns(pl.col(left_col).fill_null(pl.col(right_col))).drop(right_col)

def generate_truly_unique_cluster_id(existing_ids, denylist):
    if denylist is not None:
        with open(denylist, "r", encoding="utf-8") as f:
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
        with open(amatrix, "r", encoding="utf-8") as distance_matrix:
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
        with open(bmatrix, "r", encoding="utf-8") as distance_matrix:
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
