VERSION = "0.4.0" # does not necessarily match Tree Nine git version
print(f"PROCESS CLUSTERS - VERSION {VERSION}")

# pylint: disable=too-many-statements,too-many-branches,simplifiable-if-expression,too-many-locals,too-complex,consider-using-tuple,broad-exception-caught
# pylint: disable=wrong-import-position,useless-suppression,multiple-statements,line-too-long,consider-using-sys-exit,duplicate-code
#
# Notes:
# * Eventually we may want persistent_cluster_meta to contain parent-child cluster IDs because that might
#   prevent hypothetical edge cases where a cluster itself doesn't change but its subcluster gets renamed,
#   resulting in the parent cluster not linking to the correct MR project anymore?
# * This script calls a persistent cluster script written by Marc Perry, which handles all the tricky logic for
#   assigning persistent cluster IDs to clusters that already exist. However, we also need to assign IDs to new
#   clusters, link parent-child clusters, and upload to Microreact, which is what all this Python does.
# * There may be edge cases where Marc's script's assignment of persistent cluster IDs is non-deterministic
# * My script's assingment of brand-new cluster IDs is likely non-deterministic as it relies on sets and
#   unsorted polars dataframes. Additionally, if typical methods for assigning cluster IDs fail due to name
#   conflicts, my script will start calling random numbers to generate new cluster IDs.
# * This script is not super optimized, but it is performant (~1 minute) on laptops up to at least 4000 clusters
# * Some versions of polars are stricter than others in reading/writing TSVs and JSONs
# * 5 SNP clusters always have a 10 SNP parent, and 10 SNP clusters always have a 20 SNP parent
# * Persistent clusters can run into a Ship of Theseus situation over time 
#
# Surprisingly important information r/e logging:
#   As of mid-2025, after changing its backend to GCP Batch, Terra seems to have an issue where logging slows
# task execution to an extreme degree. We're talking "a task used to take 10 hours but now is less than
# halfway through at 48 hours on comparable inputs" levels of slowness. It genuinely appears to be faster to
# save intermediate dataframes as JSONs to the disk rather than print them to stdout, so that's what this script
# does via debug_logging_handler_df().
#   Tradeoff: If this script crashes (or if Terra/GCP has an intermittent error that crashes the VM, which seems
# to happen about 5% of the time as of mid-2025), Terra may or may not delocalize these JSON files, potentially
# leaving you without valuable debug information. If this happens, try using whatever files you have from
# find_clusters.py and/or the previous task to do testing locally.

import io
import os
import csv
import json
import random
import logging
import argparse
from datetime import datetime, timezone
import subprocess
import requests
import polars as pl # this is overkill and takes forever to import; too bad!
from polars.testing import assert_series_equal
pl.Config.set_tbl_rows(-1)
pl.Config.set_tbl_cols(-1)
pl.Config.set_tbl_width_chars(200)
pl.Config.set_fmt_str_lengths(500)
pl.Config.set_fmt_table_cell_list_len(500)
today = datetime.now(timezone.utc) # I don't care if this runs past midnight, give everything the same day!
print(f"It's {today} in Thurles right now. Up Tipp!")
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("requests").setLevel(logging.WARNING)
max_random_id_attempts = 500 # maximum attempts to fix invalid cluster IDs

if os.path.isfile("/scripts/marcs_incredible_script_update.pl"):
    script_path = "/scripts"
elif os.path.isfile("./scripts/marcs_incredible_script_update.pl"):
    script_path = "./scripts"
else:
    raise FileNotFoundError

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
    parser.add_argument('-as', '--allsamples', type=str, required=False, help='comma-delimited list of samples to consider for clustering (if absent, will do entire tree)')
    parser.add_argument('-ls', '--latestsamples', type=str, help='TSV: latest sample information (as identified by find_clusters.py)')
    #parser.add_argument('-sm', '--samplemeta', type=str, help='TSV: sample metadata pulled from terra (including myco outs), one line per sample')
    parser.add_argument('-pcm', '--persistentclustermeta', type=str, help='TSV: persistent cluster metadata from last full run of TB-D')
    parser.add_argument('-pid', '--persistentids', type=str, help='TSV: persistent IDs from last full run of TB-D')
    parser.add_argument('-mat', '--mat_tree', type=str, help='PB: tree')
    #parser.add_argument('-cs', '--contextsamples', type=int, default=0, help="[UNUSED] int: Number of context samples for cluster subtrees")
    parser.add_argument('-cd', '--combineddiff', type=str, help='diff: Maple-formatted combined diff file, needed for backmasking')
    parser.add_argument('-dl', '--denylist', type=str, required=False, help='TXT: newline delimited list of cluster IDs to never use')
    parser.add_argument('-mr', '--yes_microreact', action='store_true', help='upload clusters to MR (requires -to)')
    parser.add_argument('-d', '--today', type=str, required=True, help='ISO 8601 date, YYYY-MM-DD')
    parser.add_argument('-v', '--verbose', action='store_true', help='enable verbose logging to stdout (warning: extremely slow on Terra)')
    parser.add_argument('--disable_decimated_failsafe', action='store_true', help='do not error if a cluster on MR becomes decimated')
    parser.add_argument('--no_cleanup', action='store_true', help="do not clean up input files (this may break delocalization on Terra; only use this for rapid debug runs)")
    parser.add_argument('--mr_blank_template', type=str, help="JSON: template file for blank MR projects")
    parser.add_argument('--mr_update_template', type=str, help="JSON: template file for in-use MR projects")
    parser.add_argument('--skip_perl', action='store_true', help="skip the perl scripts to debug using existing rosetta_20/10/5 files (don't enable this for real runs!)")

    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)
    if args.persistentclustermeta and not args.persistentids:
        raise ValueError("You provided --persistentclustermeta but no --persistentids, you need both or neither")
    if args.persistentids and not args.persistentclustermeta:
        raise ValueError("You provided --persistentids but no --persistentclustermeta, you need both or neither")
    if not args.persistentids and not args.persistentclustermeta:
        start_over = True
        print("You have not provided persistent IDs nor persistent cluster metadata. This will restart clustering.")
    else:
        start_over = False
    if args.yes_microreact and not (args.mr_blank_template and args.mr_update_template):
        raise ValueError("You said --yes_microreact but didn't include --mr_blank_template and/or --mr_update_template")

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
    debug_logging_handler_df("Loaded all_latest_samples", all_latest_samples, "1_inputs")
    if not start_over:
        debug_logging_handler_df("Loaded all_persistent_samples", all_persistent_samples, "1_inputs")

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
        # TODO: this method of detecting decimation isn't working correctly and should probably be replaced
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
        debug_logging_handler_txt("Set of all latest samples", "1_inputs", 10)
        debug_logging_handler_txt(all_latest_samples_set, "1_inputs", 10)
        debug_logging_handler_txt("Set of all persistent samples", "1_inputs", 10)
        debug_logging_handler_txt(all_persistent_samples_set, "1_inputs", 10)
        if all_persistent_samples_set.issubset(all_latest_samples_set):
            debug_logging_handler_txt("All persistent samples is a subset of all latest samples", "1_inputs", 20)
        else:
            samples_missing_from_latest = all_persistent_samples_set - all_latest_samples_set # these are sets so this excludes samples exclusive to all_latest
            debug_logging_handler_txt(f"Samples appear to be missing from the latest run: {samples_missing_from_latest}", "1_inputs", 30)
            if args.allsamples:
                all_input_samples_including_unclustered = args.allsamples.split(',')
            else:
                all_input_samples_including_unclustered = None
                debug_logging_handler_txt("Missing args.allsamples; can't be sure if missing samples are dropped because they no longer cluster or if they were never input.", "1_inputs", 30)
            for sample in samples_missing_from_latest:
                if all_input_samples_including_unclustered is None:
                    pass
                elif sample in all_input_samples_including_unclustered:
                    debug_logging_handler_txt(f"{sample} is newly unclustered", "1_inputs", 30)
                else:
                    debug_logging_handler_txt(f"{sample} seems to have been dropped from inputs", "1_inputs", 30)
                # get persistent cluster ID regardless
                cluster_ids = get_cluster_ids_for_sample(all_persistent_samples, sample)
                for cluster in cluster_ids:
                    if len(get_other_samples_in_cluster(all_persistent_samples, cluster, samples_missing_from_latest)) <= 1:
                        # In theory we could handle this, in practice it's a massive pain in the neck and very easy to mess up!!
                        debug_logging_handler_txt(f"{cluster} is decimated thanks to losing all samples (or all but one)", "1_inputs", 30)
                        # IF AND ONLY IF this is not on MR (which should only happen if this is a 20-cluster with no subclusters),
                        # we can live with this being a decimated cluster.
                        if not has_microreact_url(persistent_clusters_meta, cluster):
                            debug_logging_handler_txt(f"{cluster} already lacks a Microreact URL, so we can live with it being decimated", "1_inputs", 20)
                        elif args.disable_decimated_failsafe:
                            debug_logging_handler_txt(f"{cluster} has an MR URL but we will accept it being decimated due to --disable_decimated_failsafe", "1_inputs", 30)
                        else:
                            debug_logging_handler_txt(f"{cluster} has an MR URL and should never be decimated. Cannot continue.", "1_inputs", 40)
                            exit(55)
                    else:
                        debug_logging_handler_txt(f"Dropped {sample} from {cluster} but that seems to be okay", "1_inputs", 20)
        
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
        debug_logging_handler_txt("Preparing to run the absolute legend's script...", "2_marc", 20)
        filtered_latest_20 = all_latest_20.join(all_persistent_20.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
        filtered_latest_10 = all_latest_10.join(all_persistent_10.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
        filtered_latest_5 = all_latest_5.join(all_persistent_5.drop(['cluster_id']), on="sample_id", how="inner").rename({'latest_cluster_id': 'cluster_id'}).sort('cluster_id')
        filtered_persistent_20 = all_persistent_20.join(all_latest_20.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
        filtered_persistent_10 = all_persistent_10.join(all_latest_10.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
        filtered_persistent_5 = all_persistent_5.join(all_latest_5.drop(['latest_cluster_id']), on="sample_id", how="inner").sort('cluster_id')
        for distance, dataframe in {20: filtered_latest_20, 10: filtered_latest_10, 5: filtered_latest_5}.items():
            debug_logging_handler_df(f"Filtered latest {distance}", dataframe, "2_marc")
        for distance, dataframe in {20: filtered_persistent_20, 10: filtered_persistent_10, 5: filtered_persistent_5}.items():
            debug_logging_handler_df(f"Filtered persistent {distance}", dataframe, "2_marc")

        filtered_latest_20.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_20.tsv', separator='\t', include_header=False)
        filtered_latest_10.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_10.tsv', separator='\t', include_header=False)
        filtered_latest_5.select(["sample_id", "cluster_id"]).write_csv('filtered_latest_5.tsv', separator='\t', include_header=False)
        filtered_persistent_20.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_20.tsv', separator='\t', include_header=False)
        filtered_persistent_10.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_10.tsv', separator='\t', include_header=False)
        filtered_persistent_5.select(["sample_id", "cluster_id"]).write_csv('filtered_persistent_5.tsv', separator='\t', include_header=False)

        if not args.skip_perl:
            debug_logging_handler_txt("Actually running scripts...", "2_marc", 20)
            perl_20 = subprocess.run(f"perl {script_path}/marcs_incredible_script_update.pl filtered_persistent_20.tsv filtered_latest_20.tsv", shell=True, check=True, capture_output=True, text=True)
            debug_logging_handler_txt(perl_20.stdout, "2_marc", 20)
            subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_20.tsv", shell=True, check=True)
            perl_10 = subprocess.run(f"perl {script_path}/marcs_incredible_script_update.pl filtered_persistent_10.tsv filtered_latest_10.tsv", shell=True, check=True, capture_output=True, text=True)
            debug_logging_handler_txt(perl_10.stdout, "2_marc", 20)
            subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_10.tsv", shell=True, check=True)
            perl_5 = subprocess.run(f"perl {script_path}/marcs_incredible_script_update.pl filtered_persistent_5.tsv filtered_latest_5.tsv", shell=True, check=True, capture_output=True, text=True)
            debug_logging_handler_txt(perl_5.stdout, "2_marc", 20)
            subprocess.run("mv mapped_persistent_cluster_ids_to_new_cluster_ids.tsv rosetta_stone_5.tsv", shell=True, check=True)

        # TODO: why are were we not running equalize tabs except when logging is debug?

        # debug print basic rosetta stones
        if logging.root.level == logging.DEBUG:
            for rock in ['rosetta_stone_20.tsv', 'rosetta_stone_10.tsv', 'rosetta_stone_5.tsv']:
                with open(rock, 'r', encoding="utf-8") as file:
                    debug_logging_handler_txt(f"---------------------\nContents of {rock} (before strip_tsv and equalize_tabs):\n", "2_marc", 10)
                    debug_logging_handler_txt(list(file), "2_marc", 10)
                    #subprocess.run(f"/bin/bash {script_path}/equalize_tabs.sh {rock}", shell=True, check=True)
                    
        # get more information about merges... if we have any!
        rock_pairs = {'rosetta_stone_20.tsv':'rosetta_stone_20_merges.tsv', 
                    'rosetta_stone_10.tsv':'rosetta_stone_10_merges.tsv', 
                    'rosetta_stone_5.tsv':'rosetta_stone_5_merges.tsv'}
        for rock, merge_rock in rock_pairs.items():
            if os.path.isfile(merge_rock):
                debug_logging_handler_txt(f"Found {merge_rock}, indicating clusters merged at this distance", "2_marc", 20)
                debug_logging_handler_txt(f"---------------------\nContents of {merge_rock} (before strip_tsv and equalize_tabs):\n", "2_marc", 10)
                with open(merge_rock, 'r', encoding="utf-8") as file:
                    debug_logging_handler_txt(list(merge_rock), "2_marc", 10)
                    #subprocess.run(f"/bin/bash {script_path}/equalize_tabs.sh {rock}", shell=True, check=True)
                    subprocess.run(f"/bin/bash {script_path}/strip_tsv.sh {rock} {merge_rock}", shell=True, check=True)
            else:
                debug_logging_handler_txt(f"Did not find {merge_rock}, indicating clusters didn't merge at this distance", "2_marc", 20)

        # we need schema_overrides or else cluster IDs can become non-zfilled i64
        # For some godforesaken reason, some versions of polars will throw `polars.exceptions.ComputeError: found more fields than defined in 'Schema'` even if we set
        # infer_schema = True with a hella large infer_schema_length. Idk why because the exact same file works perfectly fine on my local installation of polars (polars==1.27.0)
        # without even needing to set anything with infer_schema!! Not even a try-except with the except having a three column schema works!! Ugh!!!
        # TODO: is this because the docker is polars==1.26.0?
        # ---> WORKAROUND: equalize_tabs.sh
        debug_logging_handler_txt("Processing perl outputs...", "2_marc", 20)
        rosetta_20 = pl.read_csv("rosetta_stone_20.tsv", separator="\t", has_header=False,
            schema_overrides={"column_1": pl.Utf8, "column_2": pl.Utf8, "column_3": pl.Utf8}, 
            truncate_ragged_lines=True, ignore_errors=True, infer_schema_length=5000).rename(
            {'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id', 'column_3': 'special_handling'}
        )
        rosetta_10 = pl.read_csv("rosetta_stone_10.tsv", separator="\t", has_header=False,
            schema_overrides={"column_1": pl.Utf8, "column_2": pl.Utf8, "column_3": pl.Utf8}, 
            truncate_ragged_lines=True, ignore_errors=True, infer_schema_length=5000).rename(
            {'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id', 'column_3': 'special_handling'}
        )
        rosetta_5 = pl.read_csv("rosetta_stone_5.tsv", separator="\t", has_header=False, 
            schema_overrides={"column_1": pl.Utf8, "column_2": pl.Utf8, "column_3": pl.Utf8}, 
            truncate_ragged_lines=True, ignore_errors=True, infer_schema_length=5000).rename(
            {'column_1': 'persistent_cluster_id', 'column_2': 'latest_cluster_id', 'column_3': 'special_handling'}
        )

        # It seems theoretically possible that a (say) 20 SNP cluster could generate a persistent ID that matches a persistent ID
        # already being used by (say) 10 SNP cluster. We'll call this "cross-distance ID sharing" because I love naming things.
        # This function takes in persistent_clusters_meta so it can account for decimated cluster IDs (in theory).
        # TODO: there should be a check like this in the ad-hoc case too, just in case find_clusters does an oopsies

        debug_logging_handler_txt("Checking for cross-distance ID shares", "2_marc", 20)
        rosetta_20, rosetta_10, rosetta_5 = fix_cross_distance_ID_shares(rosetta_20, rosetta_10, rosetta_5, persistent_clusters_meta, "marc_perry")

        # TODO: Because we merge on latest_cluster_id here, and we only fixed the persistent ID, this merge could get funky?
        # In theory everything should be fine...

        debug_logging_handler_txt("Joining all_latest_samples on rosetta_20 upon latest_cluster_id to generate persistent_20_cluster_id column...", "2_marc", 20)
        latest_samples_translated = (all_latest_samples.join(rosetta_20, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_20_cluster_id'}).drop("latest_cluster_id_right")
        debug_logging_handler_txt("Joining all_latest_samples on rosetta_10 upon latest_cluster_id to generate persistent_10_cluster_id column...", "2_marc", 20)
        latest_samples_translated = (latest_samples_translated.join(rosetta_10, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_10_cluster_id'}).drop("latest_cluster_id_right")
        debug_logging_handler_txt("Nullfilling the special_handling column...", "2_marc", 20)
        latest_samples_translated = nullfill_LR(latest_samples_translated, "special_handling", "special_handling_right")
        debug_logging_handler_txt("Joining all_latest_samples on rosetta_5 upon latest_cluster_id to generate persistent_5_cluster_id column...", "2_marc", 20)
        latest_samples_translated = (latest_samples_translated.join(rosetta_5, on="latest_cluster_id", how="full")).rename({'persistent_cluster_id': 'persistent_5_cluster_id'}).drop("latest_cluster_id_right")
        debug_logging_handler_txt("Nullfilling the special_handling column (again)...", "2_marc", 20)
        latest_samples_translated = nullfill_LR(latest_samples_translated, "special_handling", "special_handling_right")
        debug_logging_handler_df("Early latest_samples_translated before polars expressions", latest_samples_translated, "2_marc")
        all_latest_samples = None # stymie my silly tendency to reuse stale variables

        debug_logging_handler_txt("Marking samples that were in 20, 10, or 5 clusters previously...", "2_marc", 20)
        latest_samples_translated = latest_samples_translated.with_columns(
            pl.when(pl.col("cluster_distance") == 20)
            .then(
                pl.when(latest_samples_translated["sample_id"].is_in(all_persistent_20["sample_id"]))
                .then(True)
                .otherwise(False)
            )
            .otherwise(None)
            .alias("in_20_cluster_last_run")
        )

        latest_samples_translated = latest_samples_translated.with_columns(
            pl.when(pl.col("cluster_distance") == 10)
            .then(
                pl.when(latest_samples_translated["sample_id"].is_in(all_persistent_10["sample_id"]))
                .then(True)
                .otherwise(False)
            )
            .otherwise(None)
            .alias("in_10_cluster_last_run")
        )

        latest_samples_translated = latest_samples_translated.with_columns(
            pl.when(pl.col("cluster_distance") == 5)
            .then(
                pl.when(latest_samples_translated["sample_id"].is_in(all_persistent_5["sample_id"]))
                .then(True)
                .otherwise(False)
            )
            .otherwise(None)
            .alias("in_5_cluster_last_run")
        )

        latest_samples_translated = latest_samples_translated.with_columns(
            pl.when(latest_samples_translated["sample_id"].is_in(all_persistent_samples_set))
            .then(False)
            .otherwise(True)
            .alias("sample_brand_new")
        )

        # Previously we did 
        # ```
        # latest_samples_translated.with_columns(
        # pl.coalesce('persistent_20_cluster_id', 'persistent_10_cluster_id', 'persistent_5_cluster_id', 'latest_cluster_id')
        # .alias("cluster_id")
        # ```
        # This was done because brand new clusters are not considered in Marc's script if they don't have any samples in a cluster
        # at that distance, since we only input sample IDs that are also in the persistent list at that SNP distance. So for example,
        # if brand new sample X and old sample Y formed brand new cluster 000033 at 10 SNPs, and Y was not in a 10 SNP cluster
        # previously, Scooby-Doo would not included in the output of Marc's script. We wouldn't have a persistent ID for it, so
        # might as well just use the latest_cluster_id (aka workdir cluster ID), righ?
        #
        # But this is problematic, since 000033 could already exist as a persistent ID in use by other samples. So we'd end up with
        # something like this:
        #
        # samp | dis | workdir | cluster_id
        # ----------------------------------
        #  A   | 10  | 000015  | 000033
        #  B   | 10  | 000015  | 000033
        #  X   | 10  | 000033  | 000033
        #  Y   | 10  | 000033  | 000033
        #
        # We had a section dedicated to adjusting this, but it was confusing and I'm not fully confident it was error-proof, so I've
        # decided to make handling the "brand new cluster" (more correctly "no persistent ID") situation more explict.
        #
        latest_samples_translated = (
            latest_samples_translated.with_columns(
                pl.coalesce(pl.col(['persistent_20_cluster_id', 'persistent_10_cluster_id', 'persistent_5_cluster_id']), pl.lit("NO_PERSIS_ID"))
                .alias("cluster_id")
            ).drop(['persistent_20_cluster_id', 'persistent_10_cluster_id', 'persistent_5_cluster_id'])
        ).rename({'latest_cluster_id': 'workdir_cluster_id'})

        # Right now samples can get "renamed" in special_handling even if their cluster didn't get renamed. Let's fix that.
        latest_samples_translated = latest_samples_translated.with_columns(
            pl.when(pl.col('workdir_cluster_id') == pl.col('cluster_id'))
            .then(pl.lit('none'))  # we don't say "unchanged" since the cluster's contents may have changed, nor do we use literal None
            .otherwise(pl.col('special_handling'))
            .alias('special_handling')
        )

        print("################# (3) SPECIAL HANDLING (of new clusters) #################")
        # This section is for handling the brand-new-cluster situation, since it generated without a persistent ID, but the workdir ID
        # it generated with could overlap with an existing persistent ID. In older versions we coalsced workdir cluster ID into (persistent)
        # cluster ID in the previous section, then in this section, detected issues by checking how many workdir cluster IDs a given
        # (persistent) cluster ID had. But it was kind of cringe so now we're handling this differently.
        debug_logging_handler_txt("Handling clusters without a persistent ID (if any)", "3_new_clusters", 20)
        latest_samples_translated = latest_samples_translated.with_columns( # this assumes all no-persistent-ids are brand new clusters; for now tis okay
            pl.when(pl.col('cluster_id') == pl.lit("NO_PERSIS_ID"))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias('cluster_brand_new')
        )
        no_persistent_id_yet = latest_samples_translated.filter(pl.col('cluster_id') == pl.lit("NO_PERSIS_ID"))
        debug_logging_handler_df("samples with no persistent ID yet", no_persistent_id_yet, "3_new_clusters")
        workdir_ids_of_no_persistent_ids = set(no_persistent_id_yet.select('workdir_cluster_id').to_series().to_list())
        for possibly_problematic_id in workdir_ids_of_no_persistent_ids:
            # check for nonsense
            n_samps_in_full_df = len(latest_samples_translated.filter(pl.col('workdir_cluster_id') == pl.lit(possibly_problematic_id)))
            n_samps_in_filtered_df = len(no_persistent_id_yet.filter(pl.col('workdir_cluster_id') == pl.lit(possibly_problematic_id)))
            assert n_samps_in_full_df == n_samps_in_filtered_df
            
            # See if this is actually a problem -- does the workdir cluster ID overlap with a persistent ID?
            # By including IDs from persistent cluster meta, this should account for decimated samples too.
            all_current_persistent_cluster_ids = set(latest_samples_translated.select('cluster_id').cast(pl.Utf8).to_series().to_list())
            all_previous_persistent_cluster_ids = set(persistent_clusters_meta.select('cluster_id').cast(pl.Utf8).to_series().to_list())
            all_cluster_ids = all_current_persistent_cluster_ids.union(all_previous_persistent_cluster_ids)

            # Yes, there's an overlap, let's generate a new ID and call that the persistent ID
            # (idk why we need to do set([str(possibly_problematic_id)]) instead of just set(str(possibly_problematic_id)) but we do, ugh)
            if set([str(possibly_problematic_id)]) & all_cluster_ids:
                new_id = generate_new_cluster_id(all_cluster_ids, "3_new_clusters")
                debug_logging_handler_txt(f"Workdir ID of brand new cluster {possibly_problematic_id} already exists as persistent ID, will change to {new_id}", "3_new_clusters", 20)
                latest_samples_translated = latest_samples_translated.with_columns([
                    pl.when(pl.col('workdir_cluster_id') == pl.lit(possibly_problematic_id))
                    .then(pl.lit(new_id))
                    .otherwise(pl.col('cluster_id'))
                    .alias('cluster_id'),

                    pl.when(pl.col('workdir_cluster_id') == pl.lit(possibly_problematic_id))
                    .then(pl.lit('brand new (renamed)'))
                    .otherwise(pl.col('special_handling'))
                    .alias('special_handling')
                ])
            
            # No overlap, let's use the workdir cluster ID as the persistent ID
            else:
                debug_logging_handler_txt(f"Workdir ID of brand new cluster {possibly_problematic_id} doesn't exist as persistent ID", "3_new_clusters", 20)
                latest_samples_translated = latest_samples_translated.with_columns([
                    pl.when(pl.col('workdir_cluster_id') == pl.lit(possibly_problematic_id))
                    .then(pl.col('workdir_cluster_id'))
                    .otherwise(pl.col('cluster_id'))
                    .alias('cluster_id'),

                    pl.when(pl.col('workdir_cluster_id') == pl.lit(possibly_problematic_id))
                    .then(pl.lit('brand new (no conflict)'))
                    .otherwise(pl.col('special_handling'))
                    .alias('special_handling')
                ])

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
            latest_samples_translated.sort('workdir_cluster_id'), "3_new_clusters")

    # ad-hoc (no-persistent-IDs) case
    else:
        latest_samples_translated = all_latest_samples.with_columns([
            pl.col("latest_cluster_id").alias("workdir_cluster_id"),
            pl.col("latest_cluster_id").alias("cluster_id"),
            pl.lit(True).alias("cluster_brand_new"),
            pl.lit("restart").alias("special_handling"),
            pl.lit(False).alias("in_20_cluster_last_run"),
            pl.lit(False).alias("in_10_cluster_last_run"),
            pl.lit(False).alias("in_5_cluster_last_run"),
            pl.lit(True).alias("sample_brand_new")
        ])

    print("################# (4) LINK PARENTS AND CHILDREN #################")
    # Possible ways to speed this up:
    # * more native polars expressions
    # * acting on the grouped dataframe instead of latest_samples_translated
    #
    # We actually do this twice, once on latest samples and once on grouped-by-persistent. In the future,
    # we may want to instead pass in persistent parent-child information as metadata so we don't need to
    # recalculate every time...
    debug_logging_handler_txt("Preparing to link parents and children...", "4_calc_paternity", 20)
    latest_samples_translated = latest_samples_translated.sort(["cluster_distance", "cluster_id"])
    debug_logging_handler_df("latest_samples_translated at start of step 4", latest_samples_translated, "4_calc_paternity")
    if not start_over:
        all_persistent_samples = all_persistent_samples.sort(["cluster_distance", "cluster_id"])
        debug_logging_handler_df("all_persistent_samples at start of step 4", all_persistent_samples, "4_calc_paternity")
    debug_logging_handler_txt("Linking samples...", "4_calc_paternity", 20)
    sample_map_latest = build_sample_map(latest_samples_translated)
    if not start_over:
        sample_map_previous = build_sample_map(all_persistent_samples)
    
    # TODO: This works, but I feel like there's bound to be another/better/faster way to do this using polars expressions
    debug_logging_handler_txt("Iterating latest_samples_translated's rows...", "4_calc_paternity", 20)
    parental_latest = establish_parenthood(latest_samples_translated, sample_map_latest)
    debug_logging_handler_txt(f"Generated latest parenthood list (len {len(parental_latest)} values), but won't update dataframe yet", "4_calc_paternity", 20)
    if not start_over:
        parental_previous = establish_parenthood(all_persistent_samples, sample_map_previous)
        debug_logging_handler_txt(f"Generated old parenthood list (len {len(parental_previous)} values), but won't update dataframe yet", "4_calc_paternity", 20)

    # We don't actually do the updates until after the group, because dealing with an agg'd list() column in polars is a mess
    # Also, acting on the grouped dataframe should be a little faster too (even if bigO doesn't really change)

    print("################# (5) GROUP #################")
    # In this section, we're going to be grouping by persistent cluster ID in order to perform some checks,
    # and get ready to check if clusters have been updated or not (however the final determination will rely
    # on a join, which happens after this, in order to properly catch clusters that lose samples)
    debug_logging_handler_txt("Grouping by persistent cluster ID...", "5_group", 20)
    latest_samples_translated = add_col_if_not_there(latest_samples_translated, "matrix_max")
    grouped = latest_samples_translated.group_by("cluster_id").agg(
        pl.col("cluster_distance").unique(),
        pl.col("matrix_max").unique(),
        pl.col("cluster_brand_new").unique(),
        pl.col("sample_brand_new").unique(),
        pl.col("special_handling").unique(),
        pl.col("workdir_cluster_id").unique(),
        pl.col("in_20_cluster_last_run").unique(),
        pl.col("in_10_cluster_last_run").unique(),
        pl.col("in_5_cluster_last_run").unique(),
        pl.col("sample_id").unique(),
        pl.col("sample_id").n_unique().alias("n_samples")
    )

    # Check every cluster has at least two samples (because this is based of the "latest" samples dataframe and doesn't have any
    # persistent metadata, we can do this check, since decimated clusters are excluded.)
    if not (grouped["sample_id"].list.len() >= 2).all(): 
        logging.basicConfig(level=logging.DEBUG) # effectively overrides global verbose
        debug_logging_handler_txt("Found cluster with less than two samples (decimated clusters are excluded in this check)", "5_group", 40)
        debug_logging_handler_df("ERROR clusters with less than two samples", grouped.filter(pl.col("sample_id").list.len() > 1), "5_group")
        raise ValueError('Found cluster with less than two samples (decimated clusters are excluded in this check')
    debug_logging_handler_txt("Asserted all clusters have at least two samples (this check happens before we have any info about decimated clusters)", "5_group", 20)

    # Check every cluster ID only has one workdir cluster ID (this is a relic of =<0.4.4's handling of brand new clusters and should never fire)
    if not (grouped["workdir_cluster_id"].list.len() <= 1).all(): 
        logging.basicConfig(level=logging.DEBUG) # effectively overrides global verbose
        debug_logging_handler_txt('Found non-zero number of "persistent" cluster IDs associated with multiple different workdir cluster IDs', "5_group", 40)
        debug_logging_handler_df("ERROR clusters with more than one workdir ID", grouped.filter(pl.col("workdir_cluster_id").list.len() > 1), "5_group")
        raise ValueError('Found non-zero number of "persistent" cluster IDs associated with multiple different workdir cluster IDs')
    debug_logging_handler_txt("Asserted all persistent cluster IDs only associated with one or zero workdir IDs", "5_group", 20)
    
    # Check only one distance per cluster ID (double checking cross-distance ID shares, this also should never fire)
    if not (grouped["cluster_distance"].list.len() == 1).all():
        debug_logging_handler_txt('Found non-zero number of "persistent" cluster IDs associated with multiple SNP distances', "5_group", 40)
        debug_logging_handler_df("ERROR clusters with not one distance", grouped.filter(pl.col("cluster_distance").list.len() != 1), "5_group")
        raise ValueError("Some clusters have multiple unique cluster_distance values.")
    debug_logging_handler_txt("Asserted all cluster_distance lists have a len of precisely 1", "5_group", 20)

    # Check only one type of special handling per cluster ID (in theory that would actually be okay but given how we do it it shouldn't happen)
    if not (grouped["special_handling"].list.len() == 1).all():
        debug_logging_handler_txt("Found different types of special handling in some clusters", "5_group", 40)
        debug_logging_handler_df("ERROR clusters with not one special_handling", grouped.filter(pl.col("special_handling").list.len() != 1), "5_group")
        raise ValueError("Some clusters have multiple unique special_handling values.")
    debug_logging_handler_txt("Asserted all special_handling lists have a len of precisely 1", "5_group", 20)

    # Check... you get the picture
    if not (grouped["cluster_brand_new"].list.len() == 1).all():
        debug_logging_handler_txt("Found clusters that don't know if they're new or not", "5_group", 40)
        debug_logging_handler_df("ERROR clusters with not one cluster_brand_new", grouped.filter(pl.col("cluster_brand_new").list.len() != 1), "5_group")
        raise ValueError("Some clusters have multiple unique cluster_brand_new values.")
    debug_logging_handler_txt("Asserted all cluster_brand_new lists have a len of precisely 1", "5_group", 20)

    debug_logging_handler_txt("Converting lists to base types where possible...", "5_group", 20)
    grouped = grouped.with_columns([
        pl.col("workdir_cluster_id").list.get(0).alias("workdir_cluster_id"),
        pl.col("cluster_distance").list.get(0).alias("cluster_distance"),
        pl.col("special_handling").list.get(0).alias("special_handling"),
        pl.col("cluster_brand_new").list.get(0).alias("cluster_brand_new")
    ])

    # Collapse the 20/10/5 last run columns
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
    ).sort('cluster_id').drop(['in_20_cluster_last_run', 'in_10_cluster_last_run', 'in_5_cluster_last_run']) # will be readded upon join

    debug_logging_handler_df("After grouping and then intager-a-fy", grouped, "5_group")
    # Previously we dropped "sample_id" column here since we grouped a second time before joining on the persistent metadata/groupby files,
    # but there isn't a reason to do that anymore.

    print("################# (6) UPDATE PATERNITY #################")
    # We already identified parents and children earlier, but now we're going to actually update the dataframe with the "updates" lists
    debug_logging_handler_txt("Updating latest grouped dataframe with paternity information...", "6_update_paternity", 20)
    grouped = grouped.with_columns(
        pl.lit(None).cast(pl.Utf8).alias("cluster_parent"),
        pl.lit([]).cast(pl.List(pl.Utf8)).alias("cluster_children")
    ).sort(["cluster_distance", "cluster_id"])
    for cluster_id, col, value in parental_latest:
        #debug_logging_handler_txt(f"For cluster {cluster_id}, col {col}, val {value} in updates", "6_update_paternity", 10) # too verbose even for debug logging
        if col == "cluster_parent":
            # Beware: This is an overwrite, so we can't check if there's multiple cluster parents
            grouped = update_cluster_column(grouped, cluster_id, "cluster_parent", value)
        else:
            grouped = grouped.with_columns(
                pl.when(pl.col("cluster_id") == cluster_id)
                .then((pl.col("cluster_children").list.concat(pl.lit(value))).list.unique())
                .otherwise(pl.col("cluster_children"))
                .alias("cluster_children")
            )
    cluster_id, parental_latest = None, None
    debug_logging_handler_df("grouped after linking parents and children", grouped, "6_update_paternity")
    if not start_over:
        debug_logging_handler_txt("Updating previous run's dataframe with paternity information...", "6_update_paternity", 20)
        persis_groupby_cluster = persis_groupby_cluster.with_columns(
            pl.lit(None).cast(pl.Utf8).alias("cluster_parent"),
            pl.lit([]).cast(pl.List(pl.Utf8)).alias("cluster_children")
        ).sort(["cluster_distance", "cluster_id"])
        for cluster_id, col, value in parental_previous:
            #debug_logging_handler_txt(f"For cluster {cluster_id}, col {col}, val {value} in updates", "6_update_paternity", 10) # too verbose even for debug logging
            if col == "cluster_parent":
                # Beware: This is an overwrite, so we can't check if there's multiple cluster parents
                persis_groupby_cluster = update_cluster_column(persis_groupby_cluster, cluster_id, "cluster_parent", value)
            else:
                persis_groupby_cluster = persis_groupby_cluster.with_columns(
                    pl.when(pl.col("cluster_id") == cluster_id)
                    .then((pl.col("cluster_children").list.concat(pl.lit(value))).list.unique())
                    .otherwise(pl.col("cluster_children"))
                    .alias("cluster_children")
                )
        cluster_id, parental_previous = None, None
        debug_logging_handler_df("persis_groupby_cluster after linking parents and children", persis_groupby_cluster, "6_update_paternity")

    # Checks involving parent/child relationships
    # The cluster_children check might change across versions of polars; right now we expect an empty list, as opposed to pl.Null or [pl.Null].
    # Previously I'm pretty sure we had [pl.Null] since we inserted paternity before the group, and now we do it after.
    # We don't use list len() because [null] is considered to have a length of 1 in some versions of polars but perhaps not others;
    # see also https://github.com/pola-rs/polars/issues/18522
    if start_over:
        check_dfs = [grouped]
    else:
        check_dfs = [grouped, persis_groupby_cluster]
    for df in check_dfs:
        assert ((df.filter(pl.col("cluster_distance") == pl.lit(5)))["cluster_parent"].is_not_null()).all(), "5-cluster with null cluster_parent"
        assert ((df.filter(pl.col("cluster_distance") == pl.lit(10)))["cluster_parent"].is_not_null()).all(), "10-cluster with null cluster_parent"
        assert ((df.filter(pl.col("cluster_distance") == pl.lit(20)))["cluster_parent"].is_null()).all(), "20-cluster with non-null cluster_parent"
        assert ((df.filter(pl.col("cluster_distance") == pl.lit(5)))["cluster_children"] == []).all(), "5-cluster with cluster_children"
        debug_logging_handler_txt("Asserted no 5 clusters have children or no parent, no 10s lack parent, and no 20s have parent", "6_update_paternity", 20)

    # convert [null] to null
    # Actually, we don't do this anymore, because I want to use pl.col("col_a").list.unique().list.sort() after joining with the 
    # persistent grouped by dataframe, in order to see if a cluster got new/different children. I'm not confident that will
    # work problem on null, so we're returning to empty lists.
    
    # debug_logging_handler_txt("Converting [null] to null in cluster_children...", "6_update_paternity", 20)
    # grouped = grouped.with_columns([
    #     # previously: pl.when(pl.col("cluster_children").list.get(0).is_null())
    #     # We used to handle paternity before grouping, resulting in empty children being [pl.Null], but now we handle
    #     # paternity after grouping so empty children are now []. In fact, list.get(0) will error in our current version. 
    #     pl.when(pl.col("cluster_children") == pl.lit([]))
    #     .then(None)
    #     .otherwise(pl.col("cluster_children"))
    #     .alias("cluster_children")
    # ])
    # if not start_over:
    #     persis_groupby_cluster = persis_groupby_cluster.with_columns([
    #         pl.when(pl.col("cluster_children") == pl.lit([]))
    #         .then(None)
    #         .otherwise(pl.col("cluster_children"))
    #         .alias("cluster_children")
    #     ])
    

    print("################# (7) JOIN with persistent information #################")
    # First, we join with the persistent cluster metadata TSV to get first_found, last_update, jurisdictions, and microreact_url
    # Then, we join with persis_groupby_cluster (which will tell us what samples clusters previously had)
    # Only after doing these can we confidentally declare which clusters have actually been updated in some way
    #
    # TODO: eventually latest cluster metadata file should be joined here too <--- nah
    if start_over:
        debug_logging_handler_txt("Generating metadata fresh (since we're starting over)...", "7_join", 20)
        all_cluster_information = grouped.with_columns([
            pl.lit(today.isoformat()).alias("first_found"), 
            pl.lit(True).alias("needs_updating"),
            pl.lit(today.isoformat()).alias("last_json_update")])
        debug_logging_handler_df("after adding relevant information", all_cluster_information, "7_join")
    else:
        debug_logging_handler_txt("Joining with the persistent metadata TSV...", "7_join", 20)
        persistent_clusters_meta = persistent_clusters_meta.with_columns(pl.lit(False).alias("cluster_brand_new"))
        all_cluster_information = grouped.join(persistent_clusters_meta, how="full", on="cluster_id", coalesce=True)

        # TODO: this is gonna require we filter out the Nones for new and decimated samples; probably isn't worth the hassle
        #assert_series_equal(
        #    all_cluster_information.filter().select("cluster_brand_new").to_series(), 
        #    all_cluster_information.filter().select("cluster_brand_new_right").to_series(),
        #    check_names=False, check_order=True
        #)

        # Persistent clutter meta can introduce decimated clusters which will have nulls for some columns, better
        # deal with those now
        all_cluster_information = all_cluster_information.with_columns(
            pl.col("cluster_brand_new").fill_null(False)
        )
        all_cluster_information = all_cluster_information.with_columns(
            pl.col("n_samples").fill_null(0)
        )

        # Set when a cluster was first found (if possible)
        all_cluster_information = all_cluster_information.with_columns([
            pl.when(pl.col("first_found").is_null())
            .then(
                pl.when(pl.col("cluster_brand_new"))
                .then(pl.lit(today.isoformat()))
                .otherwise(pl.lit("UNKNOWN")) # going foward this shouldn't happen but it did happen on older versions
            )
            .otherwise(pl.col("first_found"))
            .alias("first_found"),
        ])
        # Warn about stuff that has an unknown find date -- this shouldn't happen going forward but it did happen in the past
        # (which is why it's just a warning and not an error; for the time being I'm testing with the old JSONs)
        no_date = all_cluster_information.filter(pl.col("first_found") == pl.lit("UNKNOWN"))
        if len(no_date) > 0:
            debug_logging_handler_txt(f"Found {no_date.shape[0]} clusters with no clear first_found date", "7_join", 30)
            debug_logging_handler_df("WARNING no first_found date", no_date, "7_join")
        debug_logging_handler_df("after joining with persistent_clusters_meta", all_cluster_information, "7_join")

        # Now joined by the grouped-by-cluster persistent cluster ID information, which gives us the list of samples the clusters previously had
        # persis_groupby_cluster is only created if start_over is false, so we don't need to worry about unassigned vars here
        debug_logging_handler_txt("Joining persis_groupby_cluster...", "7_join", 20)
        all_cluster_information = all_cluster_information.join(persis_groupby_cluster, how="full", on="cluster_id", coalesce=True) # pylint: disable=possibly-used-before-assignment
        all_cluster_information = all_cluster_information.with_columns(
             pl.when(pl.col("cluster_distance").is_null())
            .then(
                pl.when(pl.col("cluster_distance_right").is_null())
                .then(pl.lit(None)) # should never happen, except in older decimated clusters
                .otherwise(pl.col("cluster_distance_right"))
            )
            .otherwise(pl.col("cluster_distance"))
            .alias("cluster_distance")
        ).drop("cluster_distance_right")
        all_cluster_information = all_cluster_information.rename({
            'sample_id_right': 'sample_id_previously',
            'cluster_parent_right': 'cluster_parent_previously',
            'cluster_children_right': 'cluster_children_previously'
        })

        # It's okay if cluster_parent is null, but the way we detect cluster children changes might freak out with nulls
        for column in ["cluster_children", "cluster_children_previously"]:
            all_cluster_information = all_cluster_information.with_columns(
                pl.col(column).fill_null([])
            )
        all_cluster_information = all_cluster_information.drop("cluster_brand_new_right")

        # Older cluster JSONs don't have a "decimated" column, so we're not gonna rely on it at all
        debug_logging_handler_txt("Declaring clusters decimated, or not...", "7_join", 20)
        if "decimated" in all_cluster_information.columns:
            all_cluster_information = all_cluster_information.drop("decimated")
        all_cluster_information = all_cluster_information.with_columns([

            # TODO: Consider adding a check for an empty list instead of just null -- seems to work though?
            pl.when(
                (
                    (pl.col('sample_id').is_null())
                    .or_(pl.col("sample_id_previously").is_null())
                )
                .and_(pl.col("cluster_brand_new") == pl.lit(False)) # We already changed these from null to False
            )
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("decimated"),

            # We distinguish between old and newly decimated clusters since that affects if they need updating,
            # and I don't really trust polars to compare empty/null lists properly
            pl.when((pl.col('sample_id').is_null()).and_(pl.col("sample_id_previously").is_not_null()))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("newly_decimated"),

            # This means a decimated cluster's persistent ID got reused, which should never happen
            pl.when(
                (pl.col('sample_id').is_not_null())
                .and_(pl.col("sample_id_previously").is_null())
                .and_(pl.col("cluster_brand_new") == pl.lit(False)) # We already changed these from null to False
            )
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("reused_decimated_persistent_id")
        ])
        reused_decimated = all_cluster_information.filter(pl.col("reused_decimated_persistent_id"))
        if len(reused_decimated) > 0:
            logging.basicConfig(level=logging.DEBUG) # effectively overrides global verbose in order to force dumping df to stdout
            debug_logging_handler_txt(f"We appear to have reused {reused_decimated.shape[0]} decimated cluster IDs", "7_join", 40)
            debug_logging_handler_df("reused decimated persistent IDs", reused_decimated, "7_join")
            raise ValueError
        all_cluster_information = all_cluster_information.drop("reused_decimated_persistent_id")
        debug_logging_handler_df("all_cluster_information at end", all_cluster_information, "7_join")

        # Now we finally have all the information we need to declare which clusters have ACTUALLY changed or not
        print("################# (8) RECOGNIZE (have I seen you before?) #################")
        # Previously we tried to get clever and rely on samples_previously_in_cluster via:
        # [True, False]/[False, True] --> some samples were in cluster previously     --> old cluster, needs updating
        # [False]                     --> no samples were in this cluster previously  --> new cluster, needs updating
        # [True]                      --> all samples were in this cluster previously --> old cluster, unchanged
        #
        # However, this doesn't work if an existing cluster splits into a new cluster where the new cluster only has old samples,
        # and likely missed other edge cases too. Still, we can detect many (most?) changes without comparing lists of samples directly. 
        #
        # Unfortunately, to be on the safe side, I still think it's worth comparing lists. I don't think there's ever been
        # a situation where the only way to have caught it is with a list compare, but it hypothetically possible.
        # The clearest way to do this in polars is to iterate the dataframe rowwise, extract the two lists as sets, and 
        # then do a set comparison. There might be a more effecient way to do this, but let's keep it simple for now.
        debug_logging_handler_txt("Determining which clusters are brand new and/or need updating...", "8_recognize", 20)
        all_cluster_information = add_col_if_not_there(all_cluster_information, "changes")
        # KEEP IN MIND: 
        # * These cases aren't mutually exclusive
        # * The easiest way to check if a polars list col is [False, True] or [True, False] is by checking its length and praying there's no nulls

        # Existing parent gains or loses children (actual child IDs checked, not just number)
        all_cluster_information = all_cluster_information.with_columns(
            pl.when(
                (pl.col("cluster_children").list.unique().list.sort() != pl.col("cluster_children_previously").list.unique().list.sort())
                .and_(pl.col('cluster_brand_new') == pl.lit(False))
            )
            .then(True)
            .otherwise(False)
            .alias("different_children")
        )
        different_children = all_cluster_information.filter(pl.col('different_children')).select(
            ['cluster_id', 'cluster_distance', 'sample_brand_new', 
            'cluster_children', 'cluster_children_previously', 
            'sample_id', 'sample_id_previously']
        )
        debug_logging_handler_txt(f"Found {different_children.shape[0]} clusters with different children", "8_recognize", 20)
        debug_logging_handler_df("different_children", different_children, "8_recognize")

        # Child has new parent (hypothetically possible if a 20/10 cluster splits weirdly enough)
        all_cluster_information = all_cluster_information.with_columns(
            pl.when(
                (pl.col('cluster_parent') != pl.col('cluster_parent_previously'))
                .and_(pl.col('cluster_brand_new') == pl.lit(False))
            )
            .then(True)
            .otherwise(False)
            .alias("new_parent")
        )
        new_parent = all_cluster_information.filter(pl.col('new_parent'))
        debug_logging_handler_txt(f"Found {new_parent.shape[0]} clusters with a new parent cluster", "8_recognize", 20)
        debug_logging_handler_df("new_parent", new_parent, "8_recognize")
        
        # The cluster is newly decimated (we print all decimated clusters but will only flag newly decimated as needs_updating)
        decimated = all_cluster_information.filter(pl.col("decimated"))
        new_decimated = all_cluster_information.filter(pl.col("newly_decimated"))
        debug_logging_handler_txt(f"Found {decimated.shape[0]} decimated clusters of which {new_decimated.shape[0]} are newly decimated", "8_recognize", 30)
        debug_logging_handler_df("decimated clusters", decimated, "8_recognize")

        # Existing cluster has brand new samples
        all_cluster_information = all_cluster_information.with_columns(
            pl.when(
                (pl.col('samples_previously_in_cluster').list.len() == pl.lit(2))
                .and_(pl.col('sample_brand_new').list.len() == pl.lit(2))
                .and_(
                    (pl.col('special_handling') == pl.lit("none"))
                    .or_(pl.col('special_handling') == pl.lit("renamed"))
                )
                .and_(pl.col('cluster_brand_new') == pl.lit(False))
            )
            .then(True)
            .otherwise(False)
            .alias("existing_new_samps")
        )
        existing_new_samps = all_cluster_information.filter(pl.col('existing_new_samps'))
        debug_logging_handler_txt(f"Found {existing_new_samps.shape[0]} existing clusters that got new samples", "8_recognize", 20)
        debug_logging_handler_df("existing_new_samps", existing_new_samps, "8_recognize")

        # Cluster is brand new (which may have only new samples, or old-and-new samples, or perhaps even only old samples)
        cluster_brand_new = all_cluster_information.filter(pl.col('cluster_brand_new'))
        debug_logging_handler_txt(f"Found {cluster_brand_new.shape[0]} brand new clusters", "8_recognize", 20)
        debug_logging_handler_df("cluster_brand_new", cluster_brand_new, "8_recognize")

        # Cluster's sample contents changed
        all_cluster_information = all_cluster_information.with_columns(
            pl.when(
                (pl.col("sample_id").list.unique().list.sort() != pl.col("sample_id_previously").list.unique().list.sort())
                .and_(pl.col('cluster_brand_new') == pl.lit(False))
            )
            .then(True)
            .otherwise(False)
            .alias("different_samples")
        )
        different_samples = all_cluster_information.filter(pl.col('different_samples'))
        debug_logging_handler_txt(f"Found {different_samples.shape[0]} clusters whose sample contents changed in some way", "8_recognize", 20)
        debug_logging_handler_df("different_samples", different_samples, "8_recognize")

        # These situations should all be covered by the above, but might be worth pulling out on their own eventually:
        # * The cluster itself is brand new, made up of only brand new samples (sample_brand_new = [true], special_handling = "brand new", samples_previously_in_cluster = [false])
        # * The cluster itself is brand new, made up of old and new samples
        # * The cluster itself is brand new, made up of only old samples (ample_brand_new = [false], special_handling = "brand new", samples_previously_in_cluster = [true])
        # * New sample caused previously unclustered existing sample to become a cluster (sample_brand_new = [false, true], special_handling = "brand new", samples_previously_in_cluster = [false])
        # * An existing cluster changed in some way (split, merge, split-n-merge) without taking in any brand new samples, but it got a new cluster ID
        # * The cluster previously exists, and has no new samples, but its name changed (if that's even possible)
        # * The cluster previously exists, and has no new samples, but it has different (don't just count number!) samples compared to previously

        # Handling for legacy JSONs
        if "last_update" in all_cluster_information.columns:
            fallback_update_col = "last_update"
            all_cluster_information = add_col_if_not_there(all_cluster_information, "last_json_update")
        else:
            fallback_update_col = "last_json_update"

        # FINALLY
        all_cluster_information = all_cluster_information.with_columns(
            pl.when(
                pl.col("different_children")
                .or_(pl.col("new_parent"))
                .or_(pl.col("newly_decimated"))
                .or_(pl.col("existing_new_samps"))
                .or_(pl.col("cluster_brand_new"))
                .or_(pl.col("different_samples"))
            )
            .then(True)
            .otherwise(False)
            .alias("needs_updating")
        )
        all_cluster_information = all_cluster_information.with_columns(
            pl.when("needs_updating")
            .then(pl.lit(today.isoformat()))
            .otherwise(pl.col(fallback_update_col))
            .alias("last_json_update")
        )

        if fallback_update_col == "last_update":
            all_cluster_information = all_cluster_information.drop("last_update")

    print("################# (9) GET NWK'D #################")
    # Pretty simple, but let's give it its own section for emphasis
    # Add some empty columns in the ad-hoc case -- parent_url and child_url will get added later
    all_cluster_information = add_col_if_not_there(all_cluster_information, "last_MR_update")
    all_cluster_information = add_col_if_not_there(all_cluster_information, "first_found")
    all_cluster_information = add_col_if_not_there(all_cluster_information, "jurisdictions")
    all_cluster_information = add_col_if_not_there(all_cluster_information, "sample_id_previously")
    all_cluster_information = add_col_if_not_there(all_cluster_information, "microreact_url")
    all_cluster_information = get_nwk_and_matrix_plus_local_mask(all_cluster_information, args.combineddiff).sort("cluster_id")
    debug_logging_handler_df("after getting nwk, matrix, and mask", all_cluster_information, "9_nwk")

    # hella_redundant is used for persistent IDs later... but maybe we should just replace it with an exploded version?
    hella_redundant = (latest_samples_translated.drop("cluster_distance")).join(grouped, on="cluster_id")
    debug_logging_handler_txt("Joined grouped with latest_samples_translated upon cluster_id to form hella_redundant", "9_nwk", 10)
    assert_series_equal(hella_redundant.select("workdir_cluster_id").to_series(), hella_redundant.select("workdir_cluster_id_right").to_series(), check_names=False)
    debug_logging_handler_txt("Asserted hella_redundant's workdir_cluster_id == hella_redundant's workdir_cluster_id_right", "9_nwk", 10)
    hella_redundant = hella_redundant.drop("workdir_cluster_id_right")
    grouped, latest_samples_translated = None, None
    debug_logging_handler_txt("Dropped workdir_cluster_id_right from hella_redundant, cleared grouped variable, cleared latest_samples_translated variable", "9_nwk", 10)

    hella_redundant = hella_redundant.with_columns([
        # this particular row is PER SAMPLE, not PER ClUSTER
        ~pl.coalesce(["in_20_cluster_last_run", "in_10_cluster_last_run", "in_5_cluster_last_run"]).alias("sample_newly_clustered")
    ]).drop("samples_previously_in_cluster")
    hella_redundant = hella_redundant.drop(["in_20_cluster_last_run", "in_10_cluster_last_run", "in_5_cluster_last_run"]).sort("cluster_id")

    # now let's get information as to which samples are new or old so we can highlight them
    debug_logging_handler_df("after processing what clusters and samples are brand new sorted by cluster_id", hella_redundant, "9_nwk")
    sample_level_information = hella_redundant.select(["sample_id", "cluster_distance", "cluster_id", "cluster_brand_new", "sample_newly_clustered", "sample_brand_new"])
    sample_level_information.write_csv(f'all_samples{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote all_samples{today.isoformat()}.tsv from some of hella_redundant's columns", "9_nwk", 20)
    sample_level_information.filter(pl.col('sample_brand_new')).write_csv(f'new_samples{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote new_samples{today.isoformat()}.tsv which should only have the brand new samples in it", "9_nwk", 20)
    sample_level_information = None


    # OLD CHECKS
    # assert (second_group["needs_updating"].list.len() == 1).all(), "Cluster not sure if it needs updating"
    # debug_logging_handler_txt("Asserted all len(needs_updating) == 1, but this may not catch all edge cases involving cluster updating", "7_secondgroup", 10)
    # second_group = second_group.with_columns(pl.col("needs_updating").list.get(0).alias("needs_updating_bool"))
    # second_group = second_group.drop("needs_updating").rename({"needs_updating_bool": "needs_updating"})
    # # TODO: why the hell was this ⬇️ a thing? was this to check it was NOT equal? ...should we maybe readd that actually?
    # # apparently i half-realized the brand-new-cluster problem at some point but didn't fix it properly
    # #assert_series_equal(second_group.select("cluster_brand_new").to_series(), second_group.select("has_new_samples").to_series(), check_names=False)
    # second_group.select(["cluster_id", "cluster_distance", "has_new_samples"]).write_csv(f"clusters_with_new_samples{today.isoformat()}.tsv", separator='\t')
    # second_group = second_group.drop("has_new_samples")

    #debug_logging_handler_df("after grouping hella_redundant by cluster_id, converting [null] to null, and other checks", second_group, "7_secondgroup")

    # This is needed to handle the no-persistent-IDs/start over situation gracefully 
    # notes: parent_url/child_urls get added later, "a_tree" etc was added by get_nwk_and_matrix_plus_local_mask()

    # okay, everything looks good so far. let's get some URLs!!
    # we already asserted that token is defined with yes_microreact hence possibly-used-before-assignment can be turned off there
    print("################# (10) MICROREACT #################")
    if args.yes_microreact:
        debug_logging_handler_txt("Assigning self-URLs...", "10_microreact", 20)
        for row in all_cluster_information.iter_rows(named=True):
            this_cluster_id = row["cluster_id"]
            distance = row["cluster_distance"]
            has_children = False if row["cluster_children"] is None else True
            brand_new = row["cluster_brand_new"]
            needs_updating = row["needs_updating"]
            URL = row["microreact_url"]

            # if a childless 20-cluster is brand new, it gets a found and an update date, but no MR URL
            # if a childless 20-cluster isn't new but has new samples, we change the update date (but still no MR URL)

            if brand_new:
                assert URL is None, f"{this_cluster_id} is brand new but already has a MR URL?"
                #all_cluster_information = update_first_found(all_cluster_information, this_cluster_id) # already did that earlier
                all_cluster_information = update_MR_datestamp(all_cluster_information, this_cluster_id)
                if distance == 20 and not has_children:
                    debug_logging_handler_txt(f"{this_cluster_id} is a brand-new 20-cluster with no children, not uploading", "10_microreact", 10)
                    #all_cluster_information = update_first_found(all_cluster_information, this_cluster_id) # already did that earlier
                    all_cluster_information = update_MR_datestamp(all_cluster_information, this_cluster_id)
                    all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "needs_updating", False)
                    continue
                URL = create_new_mr_project(token, this_cluster_id, args.mr_blank_template) # pylint: disable=possibly-used-before-assignment
                all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "microreact_url", URL)

            elif needs_updating:
                # later on, assert URL and first_found is not None... but for the first time don't do that!
                if URL is None:
                    if distance == 20 and not has_children:
                        debug_logging_handler_txt(f"{this_cluster_id} is an old 20-cluster with no children, not uploading", "10_microreact", 10)
                        all_cluster_information = update_MR_datestamp(all_cluster_information, this_cluster_id)
                        all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "needs_updating", False)
                        continue
                    debug_logging_handler_txt(f"{this_cluster_id} isn't brand new, but is flagged as needing an update and has no URL. Will make a new URL.", "10_microreact", 30)
                    URL = create_new_mr_project(token, this_cluster_id, args.mr_blank_template)
                    all_cluster_information = update_cluster_column(all_cluster_information, this_cluster_id, "microreact_url", URL)
                    all_cluster_information = update_MR_datestamp(all_cluster_information, this_cluster_id)
                else:
                    debug_logging_handler_txt(f"{this_cluster_id}'s URL ({URL}) seems valid, but we're not gonna check it", "10_microreact", 10)

            else:
                debug_logging_handler_txt(f"{this_cluster_id} seems unchanged", "10_microreact", 10)
                continue

        all_cluster_information = all_cluster_information.with_columns(
            pl.lit(None).alias("parent_URL"),
            pl.lit(None).alias("children_URLs") # TODO: how are gonna keep these ordered with the children column... does that even matter?
        )

        # now that everything has a URL, or doesn't need one, iterate a second time to get URLs of parents and children
        debug_logging_handler_txt("Searching for MR URLs of parents and children...", "10_microreact", 20)
        for row in all_cluster_information.iter_rows(named=True):
            this_cluster_id = row["cluster_id"]
            distance = row["cluster_distance"]
            has_children = False if row["cluster_children"] is None else True
            has_parent = False if row["cluster_parent"] is None else True
            needs_updating = row["needs_updating"]
            brand_new = row["cluster_brand_new"]
            URL = row["microreact_url"]

            # Because there is never a situation where a new child cluster pops up in a parent cluster that doesn't need to be updated,
            # and because MR URLs don't need to be updated, clusters that don't need updating don't need to know parent/child URLs.
            if not needs_updating:
                if row["sample_id"] is not None:
                    debug_logging_handler_txt(f"{this_cluster_id}@{distance} has no new samples (ergo no new children), skipping", "10_microreact", 10)
                else:
                    debug_logging_handler_txt(f"{this_cluster_id}@{distance} has no samples! This is likely a decimated cluster that lost all of its samples. We will not be updating its MR project.", "10_microreact", 30)
                    debug_logging_handler_txt("Row information of this decimated cluster:", "10_microreact", 10)
                    debug_logging_handler_txt(row, "10_microreact", 10)
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
        debug_logging_handler_df("all_cluster_information prior to upload", all_cluster_information, "10_microreact")

        # now that everything can be crosslinked, iterate one more time to actually upload
        # yes three iterations is weird, cringe even. I'm doing this to make debugging easier.
        debug_logging_handler_txt("Splitting and uploading...", "10_microreact", 10)
        for row in all_cluster_information.iter_rows(named=True):
            this_cluster_id = row["cluster_id"]
            distance = row["cluster_distance"]
            sample_id_list = row["sample_id"]
            has_children = False if row["cluster_children"] is None else True
            has_parent = False if row["cluster_parent"] is None else True
            parent_URL = row["parent_URL"] # children URLs handled differently
            cluster_parent = row["cluster_parent"] # None if !has_parent
            n_children = len(row["cluster_children"]) if has_children else -1
            needs_updating = row["needs_updating"]
            URL = row["microreact_url"]
            matrix_max = -1 if row["matrix_max"] is None else int(row["matrix_max"])
            bmatrix_max = -1 if row["b_max"] is None else int(row["b_max"])
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
                    debug_logging_handler_txt(f"{this_cluster_id}@{distance} marked as not needing updating (no new samples and/or childless 20-cluster), skipping", "10_microreact", 10)
                else:
                    debug_logging_handler_txt(f"You probably already know this, but {this_cluster_id}@{distance} has no samples!", "10_microreact", 30)
                continue
            
            with open(args.mr_update_template, "r", encoding="utf-8") as real_template_json:
                mr_document = json.load(real_template_json)

            # microreact needs all panals specified in a JSON file to be filled out, or else the workspace won't load

            # project title
            fullID = f"{first_found_shorthand}-{str(distance).zfill(2)}SNP-{this_cluster_id}"
            mr_document["meta"]["name"] = f"{fullID} updated {today.isoformat()}"
            mr_document["meta"]["description"] = f"{this_cluster_id} as automatically generated by version {VERSION} of process_clusters.py"

            # note
            markdown_note = f"### {this_cluster_id} ({distance}-SNP, {len(sample_id_list)} samples)\n*Updated {today.isoformat()}*\n\n"
            if matrix_max == 0:
                markdown_note += "*WARNING: This appears to be a tree where all branch lengths are 0. This is valid, but Microreact may not be able to render this cluster's NWK properly.*\n\n"
            elif bmatrix_max == 0:
                markdown_note += "*WARNING: Once backmasked, all branch lengths in this tree become 0. This is valid, but Microreact may not be able to render the backmasked NWK properly.*\n\n"
            markdown_note += f"First found {first_found}, UUID {this_cluster_id}, fullID {fullID}\n\n"
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
                debug_logging_handler_txt("Found metadata_combined.tsv, will use that for metadata", "10_microreact", 20)
                metadata_dict = csv.reader("./metadata_combined.tsv", delimiter="\t")
            else:
                debug_logging_handler_txt("Could not find metadata_combined.tsv, will mark as undefined per current CDPH guidelines", "10_microreact", 20)
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
            debug_logging_handler_txt(f"Metadata dictionary: {metadata_dict}", "10_microreact", 10)
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

            debug_logging_handler_txt(f"MR document for {this_cluster_id}:", "10_microreact", 10)
            debug_logging_handler_txt(f"{mr_document}", "10_microreact", 10)

            # actually upload
            assert URL is not None, f"No Microreact URL for {this_cluster_id}!" # new projects were already assigned a URL with blank template
            update_existing_mr_project(token, URL, mr_document, 0)
            if args.shareemail is not None:
                share_mr_project(token, URL, args.shareemail) 

        all_cluster_information = all_cluster_information.sort("cluster_id")
        debug_logging_handler_txt("Finished uploading to Microreact", "10_microreact", 20)
        debug_logging_handler_df("all_cluster_information after MR uploads", all_cluster_information, "10_microreact")
        new_persistent_meta = all_cluster_information.select(['cluster_id', 'first_found', 'last_json_update', 'last_MR_update', 'jurisdictions', 'microreact_url'])
    else:
        all_cluster_information = all_cluster_information.sort("cluster_id")
        debug_logging_handler_txt("Not touching Microreact. Final data table will NOT have microreact_url nor last_MR_update columns.", "10_microreact", 20)
        debug_logging_handler_df("all_cluster_information after no MR uploads", all_cluster_information, "10_microreact")
        new_persistent_meta = all_cluster_information.select(['cluster_id', 'first_found', 'last_json_update', 'jurisdictions'])

    print("################# (11) FINISHING UP #################")
    if not args.no_cleanup:
        if args.persistentclustermeta:
            os.remove(args.persistentclustermeta)
        if args.persistentids:
            os.remove(args.persistentids)
        if args.token is not None:
            os.remove(args.token)
        debug_logging_handler_txt("Deleted input persistentclustermeta, input persistentids, and input token", "11_finish", 10)
    all_cluster_information.write_ndjson(f'all_cluster_information{today.isoformat()}.json')
    debug_logging_handler_txt(f"Wrote all_cluster_information{today.isoformat()}.json", "11_finish", 20)
    
    # persistentMETA and persistentIDS are needed for subsequent runs
    new_persistent_meta.write_csv(f'persistentMETA{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote persistentMETA{today.isoformat()}.tsv", "11_finish", 20)
    new_persistent_ids = hella_redundant.select(['sample_id', 'cluster_id', 'cluster_distance'])
    new_persistent_ids.write_csv(f'persistentIDS{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote persistentIDS{today.isoformat()}.tsv", "11_finish", 20)

    # the sample \t cluster TSVs can be used to convert to annotated nextstrain format
    samp_persistent20cluster = new_persistent_ids.filter(pl.col('cluster_distance') == 20).select(['sample_id', 'cluster_id'])
    samp_persistent10cluster = new_persistent_ids.filter(pl.col('cluster_distance') == 10).select(['sample_id', 'cluster_id'])
    samp_persistent5cluster = new_persistent_ids.filter(pl.col('cluster_distance') == 5).select(['sample_id', 'cluster_id'])
    samp_persistent20cluster.write_csv(f'samp_persis20cluster{today.isoformat()}.tsv', separator='\t')
    samp_persistent10cluster.write_csv(f'samp_persis10cluster{today.isoformat()}.tsv', separator='\t')
    samp_persistent5cluster.write_csv(f'samp_persis5cluster{today.isoformat()}.tsv', separator='\t')
    debug_logging_handler_txt(f"Wrote samp_persis20cluster{today.isoformat()}.tsv, samp_persis10cluster{today.isoformat()}.tsv, and samp_persis5cluster{today.isoformat()}.tsv", "11_finish", 20)
    
    change_report = []
    debug_logging_handler_txt("Building change report...", "11_finish", 20)
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

    try:
        change_report_df = change_report_df.with_columns(pl.col("dist").cast(pl.Int32))
        change_report_df_no_twenties = change_report_df.filter(pl.col('dist') != 20)
    except pl.exceptions.InvalidOperationError: # decimated clusters with distance "ERROR!"
        change_report_df_no_twenties = change_report_df.filter(pl.col('dist') != "20")


    pl.Config.set_tbl_width_chars(200)
    with open(f"change_report_full{today.isoformat()}.txt", "w", encoding="utf-8") as full:
        with open(f"change_report_cdph{today.isoformat()}.txt", "w", encoding="utf-8") as cdph:
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
    debug_logging_handler_txt(f"Finished. Saved change report dataframe as change_report{today.isoformat()}.json", "11_finish", 20)

# Ultimately we want to have our cake and eat it too with regard to logging.
# 1) Terra might not delocalize log files if a task exits early
# 2) We want multiple log files when possible to make debugging actually feasible

def debug_logging_handler_txt(msg: str, logfile: str, loglevel=10):
    time = datetime.now(timezone.utc).strftime("%H:%M")
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
            f.write(f"[{str(time)}] {msg}\n")
    except Exception as e:
        print("Caught major logging error but will attempt to call logging.error to print relevant information before crashing")
        print(e)
        logging.error("Caught major logging error. This should never happen.")
        logging.error(f"   msg = {msg}") # pylint: disable=logging-fstring-interpolation
        logging.error(f"   logfile = {logfile}") # pylint: disable=logging-fstring-interpolation
        logging.error(f"   loglevel = {loglevel}") # pylint: disable=logging-fstring-interpolation
        logging.error("Crashing...")
        exit(1)

def debug_logging_handler_df(title: str, dataframe: pl.DataFrame, logfile: str):
    # If debug logging, dump dataframe to stdout and JSON (VERY SLOW ON TERRA)
    # If info+ logging, dump dataframe to JSON (faster on Terra but risky since may not delocalize on early exit)
    time = datetime.now(timezone.utc).strftime("%H:%M")
    json_name = logfile+"__"+title.replace(" ", "_")
    if logging.getLogger().getEffectiveLevel() == 10:
        debug_logging_handler_txt(f"Dataframe: {title} (see also JSON dump: {json_name}.json)", logfile, 10)
        logging.debug(dataframe)
    else:
        debug_logging_handler_txt(f"Dataframe: {title}", logfile, 20)
        debug_logging_handler_txt(f"  {json_name}.json", logfile, 20)
    try:
        dataframe.write_ndjson("./logs/"+json_name+'.json')
    except Exception as e: # ignore: broad-exception-caught
        # this isn't fatal because polars is just very picky sometimes  
        msg = "[%s @ %s] Logging error %s: Failed to write json version of debug dataframe, will attempt to save dataframe as text", logfile, time, e
        debug_logging_handler_txt(msg, logfile, 20)
        with open("./logs/"+logfile+".log", "a", encoding="utf-8") as f:
            f.write(title + "\n")
            f.write(dataframe)

def add_col_if_not_there(dataframe: pl.DataFrame, column: str):
    if column not in dataframe.columns:
        return dataframe.with_columns(pl.lit(None).alias(column))
    return dataframe

def build_sample_map(dataframe: pl.DataFrame):
    sample_map = {dist: {} for dist in [5, 10, 20]}
    for row in dataframe.iter_rows(named=True):
        sample_map[row["cluster_distance"]][row["sample_id"]] = row["cluster_id"]
        # for example:
        # {5: {'foo': '0003', 'bar': '0003'}, 10: {'foo': '0002', 'bar': '0002', 'bizz': '0002'}, 20: {'foo': '0001', 'bar': '0001', 'bizz': '0001'}}
        # Recall that every sample can only belong to one cluster at a given distance
    return sample_map

def establish_parenthood(dataframe: pl.DataFrame, sample_map: list):
    updates, cluster_with_known_parent = [], None
    for row in dataframe.iter_rows(named=True):
        cluster_id, one_sample, distance = row["cluster_id"], row["sample_id"], row["cluster_distance"]
        debug_logging_handler_txt(f"[{distance}] {one_sample} in cluster {cluster_id}", "4_calc_paternity", 10)
        # This is set up so we keep adding to updates whenever we find a child (even if that child isn't new),
        # but only add a cluster's parent once. This isn't perfect but it reduces dataframe changes later. This
        # is only helpful because latest_samples_translated is sorted by cluster ID.
        if distance == 5:
            parent_id = sample_map[10].get(one_sample)
            if parent_id and cluster_with_known_parent != cluster_id:
                cluster_with_known_parent = cluster_id
                updates.append((cluster_id, "cluster_parent", parent_id))
        elif distance == 10:
            parent_id = sample_map[20].get(one_sample)
            if parent_id and cluster_with_known_parent != cluster_id:
                cluster_with_known_parent = cluster_id
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
    return updates

def get_nwk_and_matrix_plus_local_mask(big_ol_dataframe: pl.DataFrame, combineddiff: str):
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "a_matrix")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "a_tree")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_matrix")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_tree")
    big_ol_dataframe = add_col_if_not_there(big_ol_dataframe, "b_max")
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
            btree = bmatrix = bmax = None
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
                        with open(f"b{this_cluster_id}.int", "r", encoding="utf-8") as bmaxfile:
                            bmax = int(bmaxfile.read().strip())
                    except subprocess.CalledProcessError as e:
                        logging.warning("[%s] Failed to generate locally-masked tree/matrix: %s", this_cluster_id, e.output)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_matrix", bmatrix)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_tree", btree)
                    big_ol_dataframe = update_cluster_column(big_ol_dataframe, this_cluster_id, "b_max", bmax)
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

def update_MR_datestamp(df: pl.DataFrame, cluster_id: str) -> pl.DataFrame:
    return df.with_columns(
        pl.when(df["cluster_id"] == cluster_id)
        .then(pl.lit(today.isoformat()))
        .otherwise(df["last_MR_update"])
        .alias("last_MR_update"))

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
            debug_logging_handler_txt(f"Updated MR project {URL}", "10_microreact", 10)
        else:
            debug_logging_handler_txt(f"Failed to update MR project {mr_url} [code {update_resp.status_code}]: {update_resp.text}", "10_microreact", 30)
            debug_logging_handler_txt("RETRYING...", "10_microreact", 20)
            update_existing_mr_project(token, mr_url, mr_document, retries)
    else:
        debug_logging_handler_txt(f"Failed to update MR project {mr_url} after multiple retries. Something's broken.", "10_microreact", 40)
        exit(1)

def share_mr_project(token, mr_url, email):
    api_url = "https://microreact.org/api/shares/add"
    params = {"id": mr_url}
    headers = {"Access-Token": token}
    data = {"emails": [email], "role": "viewer" }
    share_resp = requests.post(api_url, headers=headers, params=params, json=data, timeout=100)
    if share_resp.status_code != 200: 
        debug_logging_handler_txt(f"Failed to share MR project {mr_url} [code {share_resp.status_code}]: {share_resp.text}", "10_microreact", 40)
        debug_logging_handler_txt("NOT retrying as this is probably a permissions issue.", "10_microreact", 40)

def create_new_mr_project(token, this_cluster_id, mr_blank_template):
    with open(mr_blank_template, "r", encoding="utf-8") as temp_proj_json:
        mr_document = json.load(temp_proj_json)
    update_resp = requests.post("https://microreact.org/api/projects/create",
        headers={"Access-Token": token, "Content-Type": "application/json; charset=UTF-8"},
        params={"access": "private"},
        timeout=100,
        json=mr_document)
    if update_resp.status_code == 200:
        URL = update_resp.json()['id']
        debug_logging_handler_txt(f"{this_cluster_id} is brand new and has been assigned {URL} and a first-found date", "10_microreact", 10)
        return URL
    debug_logging_handler_txt(f"Failed to create new MR project for {this_cluster_id} [code {update_resp.status_code}]: {update_resp.text}", "10_microreact", 40)
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

# TODO: merge with generate_new_cluster_id()
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

def fix_cross_distance_ID_shares(rosetta_20: pl.DataFrame, rosetta_10: pl.DataFrame, rosetta_5: pl.DataFrame, persistent_clusters_meta: pl.DataFrame, logfile: str):
    """Recursive function that fixes bad cluster IDs one at a time"""
    cluster_ids_at_20 = set(rosetta_20.select('persistent_cluster_id').cast(pl.Utf8).to_series().to_list())
    cluster_ids_at_10 = set(rosetta_10.select('persistent_cluster_id').cast(pl.Utf8).to_series().to_list())
    cluster_ids_at_5 = set(rosetta_5.select('persistent_cluster_id').cast(pl.Utf8).to_series().to_list())
    cluster_ids_previous = set(persistent_clusters_meta.select('cluster_id').to_series().to_list()) # this includes decimated clusters
    all_cluster_ids = cluster_ids_at_20.union(cluster_ids_at_10).union(cluster_ids_at_5).union(cluster_ids_previous)

    # -------- 20 & 10 --------
    if cluster_ids_at_20 & cluster_ids_at_10:
        some_bad_id = (cluster_ids_at_20 & cluster_ids_at_10).pop()
        some_new_id = generate_new_cluster_id(all_cluster_ids, logfile)
        debug_logging_handler_txt(f"Found persistent ID {some_bad_id} at intersection of 20 and 10 IDs, will rename to {some_new_id}", logfile, 30)

        rosetta_10 = rosetta_10.with_columns([
            pl.when(pl.col('persistent_cluster_id') == some_bad_id)
              .then(pl.lit(some_new_id))
              .otherwise(pl.col('persistent_cluster_id'))
              .alias('persistent_cluster_id'),

            pl.when(pl.col('persistent_cluster_id') == some_bad_id)
              .then(pl.lit('renamed (cross-distance ID)'))
              .otherwise(pl.col('special_handling'))
              .alias('special_handling')
        ])

        rosetta_20, rosetta_10, rosetta_5 = fix_cross_distance_ID_shares(rosetta_20, rosetta_10, rosetta_5, persistent_clusters_meta, logfile)

    # -------- 20 & 5 --------
    if cluster_ids_at_20 & cluster_ids_at_5:
        some_bad_id = (cluster_ids_at_20 & cluster_ids_at_5).pop()
        some_new_id = generate_new_cluster_id(all_cluster_ids, logfile)
        debug_logging_handler_txt(f"Found persistent ID {some_bad_id} at intersection of 20 and 5 IDs, will rename to {some_new_id}", logfile, 30)

        rosetta_5 = rosetta_5.with_columns([
            pl.when(pl.col('persistent_cluster_id') == some_bad_id)
              .then(pl.lit(some_new_id))
              .otherwise(pl.col('persistent_cluster_id'))
              .alias('persistent_cluster_id'),

            pl.when(pl.col('persistent_cluster_id') == some_bad_id)
              .then(pl.lit('renamed (cross-distance ID)'))
              .otherwise(pl.col('special_handling'))
              .alias('special_handling')
        ])

        rosetta_20, rosetta_10, rosetta_5 = fix_cross_distance_ID_shares(rosetta_20, rosetta_10, rosetta_5, persistent_clusters_meta, logfile)

    # -------- 10 & 5 --------
    if cluster_ids_at_10 & cluster_ids_at_5:
        some_bad_id = (cluster_ids_at_10 & cluster_ids_at_5).pop()
        some_new_id = generate_new_cluster_id(all_cluster_ids, logfile)
        debug_logging_handler_txt(f"Found persistent ID {some_bad_id} at intersection of 10 and 5 IDs, will rename to {some_new_id}", logfile, 30)

        rosetta_5 = rosetta_5.with_columns([
            pl.when(pl.col('persistent_cluster_id') == some_bad_id)
              .then(pl.lit(some_new_id))
              .otherwise(pl.col('persistent_cluster_id'))
              .alias('persistent_cluster_id'),

            pl.when(pl.col('persistent_cluster_id') == some_bad_id)
              .then(pl.lit('renamed (cross-distance ID)'))
              .otherwise(pl.col('special_handling'))
              .alias('special_handling')
        ])

        rosetta_20, rosetta_10, rosetta_5 = fix_cross_distance_ID_shares(rosetta_20, rosetta_10, rosetta_5, persistent_clusters_meta, logfile)

    return [rosetta_20, rosetta_10, rosetta_5]


def generate_new_cluster_id(set_of_all_cluster_ids: set, logfile: str): # pylint: disable=inconsistent-return-statements
    """
    Returns a string with a new cluster ID.
    Usually, cluster IDs are zfilled strings representing ints, so find the maximum across all cluster IDs + 1 for new ID,
    but we also have a fallback that will try "000000", then it'll throw stuff at the wall to see what sticks. 
    * Checks and modifies GLOBAL max_random_id_attempts to prevent infinite loop *
    """
    global max_random_id_attempts # pylint: disable=global-statement

    try:
        numeric_ids = []
        for cid in set_of_all_cluster_ids:
            try:
                numeric_ids.append(int(cid))
            except Exception:
                pass

        if numeric_ids:
            new_int = max(numeric_ids) + 1
            new_id = str(new_int).zfill(6)
            if new_id not in set_of_all_cluster_ids:
                return new_id
        raise ValueError("No usable integer IDs")

    except ValueError:
        if "000000" not in set_of_all_cluster_ids:
            return "000000"
        if max_random_id_attempts > 0:
            max_random_id_attempts -= 1
            random_id = str(random.randint(1, 999999)).zfill(6)
            if random_id not in set_of_all_cluster_ids:
                return random_id
            debug_logging_handler_txt(f"Can't use {random_id} as a cluster ID, it's already in use...", logfile, 30)
            new_cluster_id = generate_new_cluster_id(set_of_all_cluster_ids, logfile)
            return new_cluster_id
        debug_logging_handler_txt("We can't seem to come up with a fixed cluster ID!", logfile, 40)
        debug_logging_handler_txt(f"Set of cluster IDs we must avoid: {set_of_all_cluster_ids}", logfile, 40)
        debug_logging_handler_txt("Giving up!", logfile, 40)
        exit(1)

if __name__ == "__main__":
    main()
