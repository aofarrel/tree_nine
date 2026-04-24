task find_CDPH_clusters {
	# Formally cluster_CDPH_method, now split into two tasks for easier debugging
	# find_clusters.py: Generates 20-10-5 clusters and distance matrices (normal and backmasked)
	# This might not work properly if any sample IDs contain a space
	input {
		File input_mat_with_new_samples
		String datestamp # has to be defined here for non-glob delocalization to work properly

		Boolean inteight                       = false
		Boolean only_matrix_special_samples    # arg is assumed to be passed in from Tree Nine
		File? special_samples
		
		File? persistent_ids
		File? persistent_cluster_meta
		File combined_diff_file           # used for local masking
		File? previous_run_cluster_json   # for comparisons -- currently we do this another way so this is unused

		# actually optional
		File? sample_metadata_tsv
		String? shareemail
		
		Int preempt = 0 # only set if you're doing a small test run
		Int memory = 50
		Boolean debug = true
		
		# temporary overrides
		File? override_find_clusters_script
		
	}

	Array[Int] cluster_distances = [20, 10, 5] # CHANGING THIS MIGHT BREAK THINGS!
	String arg_ieight = if inteight then "--int8" else ""
	
	command <<<
		set -eux pipefail
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting task"
		
		matUtils extract -i ~{input_mat_with_new_samples} -t A_big.nwk
		cp ~{input_mat_with_new_samples} .

		if [[ ! "~{override_find_clusters_script}" == '' ]]
		then
			rm /HOME/ash/scripts/find_clusters.py
			mv "~{override_find_clusters_script}" /HOME/ash/scripts/find_clusters.py
		fi

		CLUSTER_DISTANCES="~{sep=',' cluster_distances}"
		FIRST_DISTANCE="${CLUSTER_DISTANCES%%,*}"
		OTHER_DISTANCES="${CLUSTER_DISTANCES#*,}"
		echo "cluster distances $CLUSTER_DISTANCES"
		echo "First distance $FIRST_DISTANCE"
		echo "Other distances $OTHER_DISTANCES"

		# Turn off pipefail because find_clusters.py can return not-0 in non-error cases
		# TODO: this isn't great practice; there's better ways to handle the recursion;
		# should probably make the matrix generator its own script
		set +eo pipefail 

		# TODO: on very large runs, the size of $/samples may eventually cause issues with ARG_MAX
		# should be fine for our purposes though

		# shellcheck disable=SC2086
		if [[ "~{only_matrix_special_samples}" = "true" ]]
		then
			samples=$(< "~{special_samples}" tr -s '\n' ',' | head -c -1)
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

		mv A_big.nwk "BIGTREE~{datestamp}.nwk"
		echo "Renamed A_big.nwk to BIGTREE~{datestamp}.nwk"
		mv latest_samples.tsv "latest_samples~{datestamp}.tsv"
		echo "Renamed latest_samples.tsv to latest_samples~{datestamp}.tsv"
		mv latest_clusters.tsv "latest_clusters~{datestamp}.tsv"
		echo "Renamed latest_clusters.tsv to latest_clusters~{datestamp}.tsv"
		mv all_closest_relatives.txt "all_nearest_relatives~{datestamp}.txt"
		echo "Renamed all_closest_relatives.txt to all_nearest_relatives~{datestamp}.txt"
		mv unclustered_samples.txt "unclustered_samples~{datestamp}.txt"
		echo "Renamed all_closest_relatives.txt to unclustered_samples~{datestamp}.txt"

		if [[ -f lonely-subtree-assignments.tsv ]]
		then
			find . -maxdepth 1 \( -name "LONELY*.nwk" -o -name "lonely-subtree-assignments.tsv" -o -name "all_nearest_relatives*.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > unclustered_subtrees_etc.tar.gz
		else
			# untested hypothetical edge case
			echo "Couldn't find lonely-subtree-assignments.tsv, perhaps no samples are unclustered?"
		fi

		echo "The IDs of these clusters are random and DO NOT account for persistent cluster IDs. " > readme.txt
		echo "You'll need to run process_clusters.py to get your persistent cluster IDs!" >> readme.txt

		find . -maxdepth 1 \( -name "a*.nwk" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > acluster_trees.tar.gz
		find . -maxdepth 1 \( -name "b*.nwk" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > bcluster_trees.tar.gz
		find . -maxdepth 1 \( -name "a*_dmtrx" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > acluster_matrices.tar.gz
		find . -maxdepth 1 \( -name "b*_dmtrx" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > bcluster_matrices.tar.gz

		# shellcheck disable=SC2317
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished finding clusters"

	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev6"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		# The amount of outputs we originally had was overloading Terra, so some of these are commented out now.
		# Also, we try to avoid globbing where possible to make finding outs in Terra bucket easier since globs
		# create a folder with a randomized name, which is annoying!
		# 
		# Old method:  A = not internally masked, B = internally masked
		# Now we're switching to tarballs to avoid overloading Terra

		###### IMPORTANT FILES THAT SHOULD ALWAYS GO INTO SUBSEQUENT RUNS IF THEY EXIST ######
		File? new_samples = "new_samples" + datestamp + ".tsv"
		#File? clusterid_denylist = "clusterid_denylist.txt" # unused?

		# trees
		# there is no internally masked big tree because masking is done per-cluster
		File          bigtree_raw                           = "BIGTREE"+datestamp+".nwk"  # generated by matUtils (should match bigtree_gen)
		File          bigtree_gen                           = "a000000.nwk"               # generated by cluster script (should match bigtree_raw)
		File          cluster_subtrees_randomIDs            = "acluster_trees.tar.gz"
		File          cluster_subtrees_randomIDs_backmasked = "bcluster_trees.tar.gz"
		#Array[File]? acluster_trees = glob("a*.nwk")  # !UnnecessaryQuantifier
		#Array[File]? bcluster_trees = glob("b*.nwk")   # !UnnecessaryQuantifier

		# stuff related to unclustered samples
		File           all_nearest_relatives = "all_nearest_relatives" + datestamp + ".txt"
		File?          unclustered_samples = "unclustered_samples" + datestamp + ".txt"
		File           unclustered_subtrees_etc = unclustered_subtrees_etc.tar.gz
		#Array[File]?  unclustered_subtree_assignments = glob("*subtree-assignments.tsv")  # !UnnecessaryQuantifier
		#Array[File]?  unclustered_subtrees = glob("LONELY*.nwk")                          # !UnnecessaryQuantifier

		# distance matrices
		File?        bigtree_matrix = "a000000_dmtrx.tsv"
		File cluster_matrices_randomIDs  = "acluster_matrices.tar.gz"
		File cluster_matrices_randomIDs_backmasked = "acluster_matrices.tar.gz"
		#Array[File]? acluster_matrices = glob("a*_dmtrx.tsv")  # !UnnecessaryQuantifier
		#Array[File]? bcluster_matrices = glob("b*_dmtrx.tsv")  # !UnnecessaryQuantifier

		# general cluster stats
		Int n_big_clusters        = read_int("n_big_clusters")
		Int n_samples_in_clusters = read_int("n_samples_in_clusters")
		Int n_samples_processed   = read_int("n_samples_processed")
		Int n_unclustered         = read_int("n_unclustered")

		# debug
		File? latest_samples  = "latest_samples"+datestamp+".tsv"      # formerly intermediate_samplewise
		File? latest_clusters = "latest_clusters"+datestamp+".tsv"     # formerly intermediate_clusterwise

		# old old old
		#Array[File] abig_subtrees = glob("abig-subtree-*.nwk")
		#File? persistent_cluster_translator = "mapped_persistent_cluster_ids_to_new_cluster_ids.tsv"
		#Array[File] cluster_trees_json = glob("*.json")
		#Array[File] metadata_tsvs = glob("*.tsv")  # for auspice.us, which supports nwk
	}
}

task process_CDPH_clusters {
	# Formally cluster_CDPH_method, now split into two tasks for easier debugging and call cache savings.
	# This task runs DOWNSTREAM of find_CDPH_clusters.
	# find_clusters.py: Used to generate distance matrices
	# process_clusters.py: Persistent cluster IDs, subtrees, and MR upload
	# Any clusters that have at least one sample without a diff file will NOT be backmasked
	# This might not work properly if any sample IDs contain a space
	input {
		File input_mat_with_new_samples
		String datestamp # has to be defined here for non-glob delocalization to work properly

		# these need to also be in your template JSON, or else they won't show up in the default MR view
		String microreact_metadata_columns = "Epi_Duplication,Year_Collected,Patient_County,State,Country,Latitude,Longitude,Submitter_Facility,Submitter_Facility_Sample_ID,Sequencing_Facility"

		Boolean upload_clusters_to_microreact  = true
		Boolean disable_decimated_failsafe     = false
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
		Boolean debug = true
		
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
	String arg_disable_decimated_failsafe = if disable_decimated_failsafe then "--disable_decimated_failsafe" else ""
	
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
			echo "This will ALSO skip matrix_max calculations due to the lack of latest_clusters.tsv!"
			mv "~{override_latest_samples_tsv}" ./latest_samples.tsv
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
		echo "Args:"
		echo "--latestsamples latest_samples.tsv $LATEST_CLUSTERS_META"
		echo "-mat ~{input_mat_with_new_samples}"
		echo "-cd ~{combined_diff_file}"
		echo "--mr_metadata_columns ~{microreact_metadata_columns}"
		echo "--entity_id ~{arg_denylist} ~{arg_shareemail} ~{arg_microreact} --today ~{datestamp} ~{arg_disable_decimated_failsafe}"
		echo "--no_err_on_decimated_on_mr $TOKEN_ARG $MR_UPDATE_JSON_ARG $MR_BLANK_JSON_ARG $MR_DECIMATED_JSON_ARG"
		echo "$PERSISTENTIDS_ARG $PERSISTENTMETA_ARG $ALLSAMPLES_ARG_1 $ALLSAMPLES_ARG_2 $SAMPLEMETADATA_ARG"

		echo "Running second script"

		# shellcheck disable=SC2086 # already dquoted
		python3 /HOME/ash/scripts/process_clusters.py \
			--latestsamples latest_samples.tsv $LATEST_CLUSTERS_META \
			-mat "~{input_mat_with_new_samples}" \
			-cd "~{combined_diff_file}" \
			--mr_metadata_columns ~{microreact_metadata_columns} \
			--entity_id ~{arg_denylist} ~{arg_shareemail} ~{arg_microreact} --today ~{datestamp} ~{arg_disable_decimated_failsafe} \
			--no_err_on_decimated_on_mr $TOKEN_ARG $MR_UPDATE_JSON_ARG $MR_BLANK_JSON_ARG $MR_DECIMATED_JSON_ARG \
			$PERSISTENTIDS_ARG $PERSISTENTMETA_ARG $ALLSAMPLES_ARG_1 $ALLSAMPLES_ARG_2 $SAMPLEMETADATA_ARG 

		PY_EXIT_CODE=$? # this does not seem reliable on WDL nowadays? hmmmm...

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
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished find_clusters.py"
		fi
		if [ ~{debug} = "true" ]; then ls -lha; fi

		# if process_clusters.py errored, NOW we should crash, since we have logs and such
		exit $PY_EXIT_CODE

		# shellcheck disable=SC2317
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished task"

	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev6"
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
		# do not exist if skipping find_clusters.py
		Int n_big_clusters        = read_int("n_big_clusters")
		Int n_samples_in_clusters = read_int("n_samples_in_clusters")
		Int n_samples_processed   = read_int("n_samples_processed")
		Int n_unclustered         = read_int("n_unclustered")

		# debug
		File? logs = "logs.zip"
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