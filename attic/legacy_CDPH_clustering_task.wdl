version 1.0

# This task was split into two tasks in version 0.6.8 to enable better call-cache handling and quicker debugging.
# This all-in-one legacy version is no longer supported and should not be used.

task cluster_CDPH_method {
	# find_clusters.py: Generates 20-10-5 clusters and distance matrices (normal and backmasked)
	# process_clusters.py: Persistent cluster IDs, subtrees, and MR upload
	# Any clusters that have at least one sample without a diff file will NOT be backmasked
	# This might not work properly if any sample IDs contain a space
	input {
		File input_mat_with_new_samples
		String datestamp # has to be defined here for non-glob delocalization to work properly

		# these need to also be in your template JSON, or else they won't show up in the default MR view
		String microreact_metadata_columns = "Epi_Duplication,Year_Collected,Patient_County,State,Country,Latitude,Longitude,Submitter_Facility,Submitter_Facility_Sample_ID,Sequencing_Facility"

		Boolean upload_clusters_to_microreact  = true
		Boolean no_dropped_sample_failsafe     = false
		Boolean inteight                       = false
		Boolean force_microreact_update        = false
		Boolean only_matrix_special_samples    # arg is assumed to be passed in from Tree Nine
		File? special_samples
		
		File? persistent_denylist
		File? persistent_ids
		File? persistent_cluster_meta
		File combined_diff_file           # used for local masking
		File? previous_run_cluster_json   # for comparisons -- currently we do this another way so this is unused

		# keep these files in the workspace bucket for now
		File? microreact_decimated_template_json
		File? microreact_update_template_json
		File? microreact_blank_template_json  # hardcoded to expect a file named BLANK_template.json
		File? microreact_key

		# actually optional
		File? sample_metadata_tsv
		String? shareemail
		
		Int preempt = 0 # only set if you're doing a small test run
		Int memory = 50
		Boolean verbose = false
		
		# temporary overrides
		File? override_latest_samples_tsv  # if provided, skips find_clusters.py
		File? override_find_clusters_script
		File? override_process_clusters_script
		File? override_summarize_changes_script
		
	}
	# We cannot `String arg_token = if upload_clusters_to_microreact then "--token ~{microreact_key}" else "" ` or else the literal gs:// will
	# instead of the delocalized version, so some args will need to be handled in the command section itself

	Array[Int] cluster_distances = [20, 10, 5] # CHANGING THIS WILL BREAK SECOND SCRIPT!
	String arg_denylist = if defined(persistent_denylist) then "--dl ~{persistent_denylist}" else ""
	String arg_shareemail = if defined(shareemail) then "-s ~{shareemail}" else ""
	String arg_microreact = if upload_clusters_to_microreact then "--upload_to_microreact" else ""
	String arg_ieight = if inteight then "--int8" else ""
	String arg_disable_dropped_sample_failsafe = if no_dropped_sample_failsafe then "--no_dropped_sample_failsafe" else ""
	String arg_force_mr_update = if force_microreact_update then "--force_mr_update" else ""
	String arg_verbose = if verbose then "--verbose" else ""
	
	command <<<
	set -eux pipefail
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting task"
		
	# validate inputs
		if [[ "~{upload_clusters_to_microreact}" = "true" ]]
		then
			if [ -f "~{microreact_key}" ]
			then
				TOKEN_ARG="--token ~{microreact_key}"
			else
				echo "Upload to microreact is true, but no token provided. Crashing!"
				exit 1
			fi

			if [ -f "~{microreact_update_template_json}" ]
			then
				MR_UPDATE_JSON_ARG="--mr_update_template ~{microreact_update_template_json}"
			else
				echo "Upload to microreact is true, but no microreact_update_template_json provided. Crashing!"
				exit 1
			fi

			if [ -f "~{microreact_blank_template_json}" ]
			then
				MR_BLANK_JSON_ARG="--mr_blank_template ~{microreact_blank_template_json}"
			else
				echo "Upload to microreact is true, but no microreact_blank_template_json provided. Crashing!"
				exit 1
			fi

			if [ -f "~{microreact_decimated_template_json}" ]
			then
				MR_DECIMATED_JSON_ARG="--mr_decimated_template ~{microreact_decimated_template_json}"
			else
				echo -n "Upload to microreact is true, but no microreact_decimated_template_json provided. This isn't recommended,"
				echo -n "because decimated clusters on Microreact will never be updated, which might lead to incorrect assumptions."
			fi
		else
			TOKEN_ARG=""
			MR_UPDATE_JSON_ARG=""
			MR_BLANK_JSON_ARG=""
			MR_DECIMATED_JSON_ARG=""
		fi

		# we do similar logic within process_clusters.py too, but if we can crash before find_clusters.py that'd be ideal
		if [ -f "~{persistent_ids}" ]
		then
			if [ -f "~{persistent_cluster_meta}" ]
			then
				PERSISTENTIDS_ARG="--persistentids ~{persistent_ids}"
				PERSISTENTMETA_ARG="--persistentclustermeta ~{persistent_cluster_meta}"
			else
				echo "Found persistent IDs file but no persistent cluster meta. You need neither or both. Crashing!"
				exit 1
			fi
		else
			if [ -f "~{persistent_cluster_meta}" ]
			then
				echo "Found persistent cluster meta file but no persistent IDs. You need neither or both. Crashing!"
				exit 1
			else
				echo "Found neither persistent IDs file nor persistent cluster meta, will be running without persistent IDs"
				PERSISTENTIDS_ARG=""
				PERSISTENTMETA_ARG=""
			fi
		fi

		if [ -f "~{sample_metadata_tsv}" ]
		then
			SAMPLEMETADATA_ARG="--samplemeta ~{sample_metadata_tsv}"
		else
			SAMPLEMETADATA_ARG=""
		fi

		matUtils extract -i ~{input_mat_with_new_samples} -t A_big.nwk
		cp ~{input_mat_with_new_samples} .

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading files"

		#wget https://raw.githubusercontent.com/aofarrel/tree_nine/0.6.4/find_clusters.py
		#wget https://raw.githubusercontent.com/aofarrel/tree_nine/0.6.4/process_clusters.py
		#wget https://raw.githubusercontent.com/aofarrel/tree_nine/0.6.4/summarize_changes_alt.py

		if [[ ! "~{override_find_clusters_script}" == '' ]]
		then
			rm /HOME/ash/scripts/find_clusters.py
			mv "~{override_find_clusters_script}" /HOME/ash/scripts/find_clusters.py
		fi

		if [[ ! "~{override_process_clusters_script}" == '' ]]
		then
			rm /HOME/ash/scripts/process_clusters.py
			mv "~{override_process_clusters_script}" /HOME/ash/scripts/process_clusters.py
		fi

		if [[ ! "~{override_summarize_changes_script}" == '' ]]
		then
			rm /HOME/ash/scripts/summarize_changes_alt.py
			mv "~{override_summarize_changes_script}" /HOME/ash/scripts/summarize_changes_alt.py
		fi

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Files moved if necessary. Workdir:"
		tree

		CLUSTER_DISTANCES="~{sep=',' cluster_distances}"
		FIRST_DISTANCE="${CLUSTER_DISTANCES%%,*}"
		OTHER_DISTANCES="${CLUSTER_DISTANCES#*,}"
		echo "cluster distances $CLUSTER_DISTANCES"
		echo "First distance $FIRST_DISTANCE"
		echo "Other distances $OTHER_DISTANCES"

		# Turn off pipefail at this point for a few reasons
		# 1) find_clusters.py can return not-0 in non-error cases
		# 2) process_clusters.py writes a lot of logs to disk and we need them if it fails
		set +eo pipefail 

		# TODO: on very large runs, the size of $/samples may eventually cause issues with ARG_MAX
		# should be fine for our purposes though

		if [[ "~{only_matrix_special_samples}" = "true" ]]
		then
			samples=$(< "~{special_samples}" tr -s '\n' ',' | head -c -1)
			ALLSAMPLES_ARG_1="--allsamples"
			ALLSAMPLES_ARG_2="$samples"
		else
			ALLSAMPLES_ARG_1=""
			ALLSAMPLES_ARG_2=""
		fi

		if [[ "~{override_latest_samples_tsv}" == '' ]]
		then

			# shellcheck disable=SC2086
			if [[ "~{only_matrix_special_samples}" = "true" ]]
			then
				echo "Samples that will be in the distance matrix: $samples"
				echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running find_clusters.py"

				python3 /HOME/ash/scripts/find_clusters.py \
					"~{input_mat_with_new_samples}" \
					--samples $samples \
					--collection-name big \
					-t NB \
					-d "$FIRST_DISTANCE" \
					-rd "$OTHER_DISTANCES" \
					-v ~{arg_ieight}
			else
				echo "No sample selection file passed in, will matrix the entire tree (WARNING: THIS MAY BE VERY SLOW)"
				echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running find_clusters.py"
				python3 /HOME/ash/scripts/find_clusters.py \
					"~{input_mat_with_new_samples}" \
					--collection-name big \
					-t NB \
					-d "$FIRST_DISTANCE" \
					-rd "$OTHER_DISTANCES" \
					-v ~{arg_ieight}
			fi
			LATEST_CLUSTERS_META="--latestclustermeta latest_clusters.tsv"
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished running find_clusters.py"

		else
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipped find_clusters.py because you provided override_latest_samples_tsv"
			echo "This requires we create bogus fallback files for n_big_clusters, etc, as WDL cannot handle read_int() on non-existent file (even if var is marked optional); these will show as -1"
			echo "This will ALSO skip matrix_max calculations due to the lack of latest_clusters.tsv!"
			mv "~{override_latest_samples_tsv}" ./latest_samples.tsv
			echo "-1" > n_big_clusters
			echo "-1" > n_samples_in_clusters
			echo "-1" > n_samples_processed
			echo "-1" > n_unclustered
			LATEST_CLUSTERS_META=""
			tree

		fi

		echo "Current sample information:"
		cat latest_samples.tsv
		echo "Contents of workdir:"
		tree
		# A_big.nwk									big tree, nwk format (will be renamed later)
		# LONELY-subtree-n.nwk (n as variable)		subtrees (usually multiple) of unclustered samples
		# unclustered_samples.txt					what it says on the tin
		# lonely-subtree-assignments.tsv			which subtree each unclustered sample ended up in
		# cluster_annotation_workdirIDs.tsv			can be used to annotate by nonpersistent cluster (but isn't, at least not yet)
		# latest_samples.tsv						used by persistent ID script (will be renamed later)
		# n_big_clusters (n as constant)			# of 20SNP clusters
		# n_samples_in_clusters (n as constant)		# of samples that clustered
		# n_samples_processed (n as constant)		# of samples processed by find_clusters.py
		# n_unclustered (n as constant)				# of samples that failed to cluster
		# ...and one distance matrix per cluster, and also one(?) subtree per cluster. Later, there will be two of each per cluster thanks to backmasking

		mkdir logs

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generated these args for process_clusters.py:"
		echo "--combineddiff ~{combined_diff_file}"
		echo "--latestsamples latest_samples.tsv"
		echo "$LATEST_CLUSTERS_META"
		echo "--mat_tree ~{input_mat_with_new_samples}"
		echo "--today ~{datestamp}"
		echo "~{arg_denylist}"
		echo "~{arg_disable_dropped_sample_failsafe}"
		echo "~{arg_verbose}"
		echo "$PERSISTENTIDS_ARG $PERSISTENTMETA_ARG"
		echo "$SAMPLEMETADATA_ARG"
		
		# microreact stuff
		echo "--mr_metadata_columns ~{microreact_metadata_columns}"
		echo "~{arg_force_mr_update}"
		echo "~{arg_microreact}"
		echo "~{arg_shareemail}"
		echo "$MR_UPDATE_JSON_ARG $MR_BLANK_JSON_ARG $MR_DECIMATED_JSON_ARG"
		echo "$TOKEN_ARG"
		
		# list of samples (last because longest)
		echo "$ALLSAMPLES_ARG_1 $ALLSAMPLES_ARG_2"

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running process_clusters.py"

		# shellcheck disable=SC2086 # already dquoted
		python3 /HOME/ash/scripts/process_clusters.py \
			--combineddiff "~{combined_diff_file}" \
			--latestsamples latest_samples.tsv \
			$LATEST_CLUSTERS_META \
			--mat_tree "~{input_mat_with_new_samples}" \
			--today ~{datestamp} \
			~{arg_denylist} \
			~{arg_disable_dropped_sample_failsafe} \
			~{arg_verbose} \
			$PERSISTENTIDS_ARG $PERSISTENTMETA_ARG \
			$SAMPLEMETADATA_ARG \
			--mr_metadata_columns ~{microreact_metadata_columns} \
			~{arg_force_mr_update} \
			~{arg_microreact} \
			~{arg_shareemail} \
			$MR_UPDATE_JSON_ARG $MR_BLANK_JSON_ARG $MR_DECIMATED_JSON_ARG \
			$TOKEN_ARG \
			$ALLSAMPLES_ARG_1 $ALLSAMPLES_ARG_2 

		PY_EXIT_CODE=$? # might intermittently fails on Terra (it happened once and now I'm scared)

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Zipping process_clusters.py's logs"
		zip -r logs.zip ./logs
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Logs zipped"

		if [ -f "rosetta_stone_20_merges.tsv" ]
		then
			echo "Rosetta 20 merges"
			cat rosetta_stone_20_merges.tsv
		fi
		if [ -f "rosetta_stone_10_merges.tsv" ]
		then
			echo "Rosetta 10 merges"
			cat rosetta_stone_10_merges.tsv
		fi
		if [ -f "rosetta_stone_5_merges.tsv" ]
		then
			echo "Rosetta 5 merges"
			cat rosetta_stone_5_merges.tsv
		fi

		if [ "~{previous_run_cluster_json}" != "" ]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running summarize_changes_alt.py"
			python3 /HOME/ash/scripts/summarize_changes_alt.py "all_cluster_information~{datestamp}.json"
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished summarize_changes_alt.py"
		fi
		if [ ~{verbose} = "true" ]; then tree; fi
		

		# MR templates are generally deleted in the script itself to avoid globbing with the subtrees
		mv A_big.nwk "BIGTREE~{datestamp}.nwk"
		echo "Renamed A_big.nwk to BIGTREE~{datestamp}.nwk"
		mv latest_samples.tsv "latest_samples~{datestamp}.tsv"
		echo "Renamed latest_samples.tsv to latest_samples~{datestamp}.tsv"
		if [[ "~{override_latest_samples_tsv}" == '' ]]
		then
			mv latest_clusters.tsv "latest_clusters~{datestamp}.tsv"
			echo "Renamed latest_clusters.tsv to latest_clusters~{datestamp}.tsv"
			mv all_closest_relatives.txt "all_nearest_relatives~{datestamp}.txt"
			echo "Renamed all_closest_relatives.txt to all_nearest_relatives~{datestamp}.txt"
			mv unclustered_samples.txt "unclustered_samples~{datestamp}.txt"
			echo "Renamed all_closest_relatives.txt to unclustered_samples~{datestamp}.txt"
		else
			echo "override_latest_samples_tsv is non-null, so we won't output latest_clusters, all_nearest_relatives, or unclustered_samples"
		fi


		# if process_clusters.py errored, NOW we should crash, since we have logs and such
		exit $PY_EXIT_CODE

		# shellcheck disable=SC2317
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished task"

	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev13"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		# The amount of outputs we originally had was overloading Terra, so some of these are commented out now.
		# Also, we try to avoid globbing where possible to make finding outs in Terra bucket easier since globs
		# create a folder with a randomized name, which is annoying!

		###### IMPORTANT FILES THAT SHOULD ALWAYS GO INTO SUBSEQUENT RUNS IF THEY EXIST ######
		File? new_samples = "new_samples" + datestamp + ".tsv"
		File? clusterid_denylist = "clusterid_denylist.txt"
		File? new_persistent_ids = "persistentIDS" + datestamp + ".tsv"
		File? new_persistent_meta = "persistentMETA" + datestamp + ".tsv"
		File? final_cluster_information_json = "all_cluster_information" + datestamp + ".json"
		File? change_report_json = "change_report" + datestamp + ".json"

		File? updated_mr_URIs_file = "updated_mr_URIs" + datestamp + ".txt"

		# trees -- A = not internally masked, B = internally masked
		# there is no internally masked big tree because masking is done per-cluster
		File? bigtree_raw    = "BIGTREE"+datestamp+".nwk"   # generated directly via matUtils
		File? bigtree_gen    = "a000000.nwk"           # generated by cluster script (but should be equivalent to bigtree_raw)
		#Array[File]? acluster_trees = glob("a*.nwk")  # !UnnecessaryQuantifier
		Array[File]? bcluster_trees = glob("b*.nwk")   # !UnnecessaryQuantifier
		
		# stuff related to unclustered samples
		File?         all_nearest_relatives = "all_nearest_relatives" + datestamp + ".txt"
		Array[File]?  unclustered_subtree_assignments = glob("*subtree-assignments.tsv")  # !UnnecessaryQuantifier
		#Array[File]? unclustered_subtrees = glob("LONELY*.nwk")                          # !UnnecessaryQuantifier
		File?         unclustered_samples = "unclustered_samples" + datestamp + ".txt"

		# distance matrices
		File?        bigtree_matrix = "a000000_dmtrx.tsv"
		Array[File]? acluster_matrices = glob("a*_dmtrx.tsv")  # !UnnecessaryQuantifier
		Array[File]? bcluster_matrices = glob("b*_dmtrx.tsv")  # !UnnecessaryQuantifier

		# general cluster stats
		# marked as -1 if skipping find_clusters.py because read_int() fails if file doesn't exist (see inline echo)
		Int n_big_clusters        = read_int("n_big_clusters")
		Int n_samples_in_clusters = read_int("n_samples_in_clusters")
		Int n_samples_processed   = read_int("n_samples_processed")
		Int n_unclustered         = read_int("n_unclustered")

		# debug
		File logs = "logs.zip"
		File? change_report_full       = "change_report_full"+datestamp+".txt"  # all clusters
		File? change_report_cdph       = "change_report_cdph"+datestamp+".txt"  # excludes 20-clusters
		File? intermediate_samplewise  = "latest_samples"+datestamp+".tsv"      # from find_clusters.py
		File? intermediate_clusterwise = "latest_clusters"+datestamp+".tsv"     # from find_clusters.py, currently only for matrix_max
		File? input_metadata_tsv       = sample_metadata_tsv                    # because metadata is mutable on Terra

		# for annotation of trees
		File? samp_cluster_twn = "samp_persis20cluster" + datestamp + ".tsv" # this format is specifically for nextstrain conversion
		File? samp_cluster_ten = "samp_persis10cluster" + datestamp + ".tsv" # this format is specifically for nextstrain conversion
		File? samp_cluster_fiv = "samp_persis5cluster"  + datestamp + ".tsv" # this format is specifically for nextstrain conversion

		# old old old
		#Array[File] abig_subtrees = glob("abig-subtree-*.nwk")
		#File? persistent_cluster_translator = "mapped_persistent_cluster_ids_to_new_cluster_ids.tsv"
		#Array[File] cluster_trees_json = glob("*.json")
		#Array[File] metadata_tsvs = glob("*.tsv")  # for auspice.us, which supports nwk
	}
}
