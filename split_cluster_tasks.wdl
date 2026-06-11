version 1.0

task find_CDPH_clusters {
	# Formally cluster_CDPH_method, now split into two tasks for easier debugging
	# find_clusters.py: Generates 20-10-5 clusters and distance matrices (normal and backmasked)
	# This might not work properly if any sample IDs contain a space
	input {
		File input_mat_with_new_samples
		String datestamp # has to be defined here for non-glob delocalization to work properly

		# If not provided, it is assumed you want to matrix the ENTIRE tree
		Boolean only_matrix_special_samples
		File? special_samples

		Int memory = 50

		# EXPERIMENTAL: Store the distance matrix in memory as eight-bit unsigned integers. Actual calculations are
		# done in 64 bit and anything that would overflow is set to 255. This can resolve out-of-memory issues on
		# limited hardware, but it's not recommended.
		Boolean inteight = false
		
		# these should only be set for test runs/debugging
		File?   override_find_clusters_script
		Int     preempt = 0
		
	}

	Array[Int] cluster_distances = [20, 10, 5] # CHANGING THIS MIGHT BREAK THINGS!
	String arg_ieight = if inteight then "--int8" else ""
	
	command <<<
		set -eux pipefail
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting task"
		
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting newick (A_big.nwk) from input PB"
		matUtils extract -i ~{input_mat_with_new_samples} -t A_big.nwk
		cp ~{input_mat_with_new_samples} .

		if [[ ! "~{override_find_clusters_script}" == '' ]]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Overwriting find_clusters.py script with user-defined input"
			rm /HOME/ash/scripts/find_clusters.py
			mv "~{override_find_clusters_script}" /HOME/ash/scripts/find_clusters.py
		fi

		CLUSTER_DISTANCES="~{sep=',' cluster_distances}"
		FIRST_DISTANCE="${CLUSTER_DISTANCES%%,*}"
		OTHER_DISTANCES="${CLUSTER_DISTANCES#*,}"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cluster distances: $CLUSTER_DISTANCES"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] First distance: $FIRST_DISTANCE"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Other distances: $OTHER_DISTANCES"

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
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished running find_clusters.py"

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Contents of workdir before processing outputs:"
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
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed A_big.nwk to BIGTREE~{datestamp}.nwk"
		mv all_closest_relatives.txt "all_nearest_relatives~{datestamp}.txt"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed all_closest_relatives.txt to all_nearest_relatives~{datestamp}.txt"
		mv latest_samples.tsv "latest_samples~{datestamp}.tsv"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed latest_samples.tsv to latest_samples~{datestamp}.tsv"
		mv latest_clusters.tsv "latest_clusters~{datestamp}.tsv"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed latest_clusters.tsv to latest_clusters~{datestamp}.tsv"
		mv unclustered_samples.txt "unclustered_samples~{datestamp}.txt"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed unclustered_samples.txt to unclustered_samples~{datestamp}.txt"

		# copy stuff that will go into an archive but we also want as a task-level output
		cp aworkdir000000.nwk bigtree_gen.temp
		cp "unclustered_samples~{datestamp}.txt" nonclustered_samples.temp

		find . -maxdepth 1 \( -name "LONELY*.nwk" -o -name "lonely-subtree-assignments.tsv" -o -name "unclustered*.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > unclustered_subtrees_etc.tar.gz

		echo "The IDs of these clusters are random and DO NOT account for persistent cluster IDs. " > readme.txt
		echo "You'll need to run process_clusters.py to get your persistent cluster IDs!" >> readme.txt

		find . -maxdepth 1 \( -name "a*.nwk" -o -name "a*.pb" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > randomID_cluster_trees.tar.gz
		find . -maxdepth 1 \( -name "a*_dmtrx.tsv" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > randomID_cluster_matrices.tar.gz

		# for output matching (will include "workdir" for consistency)
		mv bigtree_gen.temp aworkdir000000.nwk
		mv nonclustered_samples.temp "unclustered_samples~{datestamp}.txt"

		# shellcheck disable=SC2317
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Contents of workdir after processing outputs:"
		tree
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished finding clusters"

	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev14"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File all_nearest_relatives    = "all_nearest_relatives" + datestamp + ".txt" # includes unclustered + clustered samples
		File latest_samples_tsv       = "latest_samples"+datestamp+".tsv"            # formerly intermediate_samplewise
		File latest_clusters_tsv      = "latest_clusters"+datestamp+".tsv"           # formerly intermediate_clusterwise
		File unclustered_samples      = "unclustered_samples" + datestamp + ".txt"
		File unclustered_subtrees_etc = "unclustered_subtrees_etc.tar.gz"            # contains subtree assignment information

		# trees and matrices
		File      bigtree_gen                           = "aworkdir000000.nwk"               # generated by cluster script (should match bigtree_raw)
		File      bigtree_matrix                        = "aworkdir000000_dmtrx.tsv"
		File      bigtree_raw                           = "BIGTREE"+datestamp+".nwk"         # generated by matUtils (should match bigtree_gen)
		File      cluster_matrices_randomIDs            = "randomID_cluster_matrices.tar.gz" # formerly Array[File]? acluster_matrices
		File      cluster_subtrees_randomIDs            = "randomID_cluster_trees.tar.gz"    # formerly Array[File]? acluster_trees
		
		# stats for the nerds
		Int n_big_clusters        = read_int("n_big_clusters")
		Int n_samples_in_clusters = read_int("n_samples_in_clusters")
		Int n_samples_processed   = read_int("n_samples_processed")
		Int n_unclustered         = read_int("n_unclustered")
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

		# These come from find_CDPH_clusters WDL task/find_clusters.py
		File latest_samples_tsv
		File latest_clusters_tsv
		File? cluster_matrices_randomIDs_tarball
		File? cluster_subtrees_randomIDs_tarball

		Array[String]? microreact_metadata_columns
		Boolean use_hardcoded_column_renames   = false

		Boolean upload_clusters_to_microreact  = true
		Boolean no_dropped_sample_failsafe     = false
		Boolean only_matrix_special_samples    # arg is assumed to be passed in from Tree Nine
		File? special_samples

		Boolean force_microreact_update        = false
		
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
		Boolean verbose = true
		
		# temporary overrides
		File? override_find_clusters_script
		File? override_process_clusters_script
		File? override_summarize_changes_script
		File? override_mass_rename_script
		
	}
	# We cannot `String arg_token = if upload_clusters_to_microreact then "--token ~{microreact_key}" else "" ` or else the literal gs:// will
	# instead of the delocalized version, so some args will need to be handled in the command section itself

	Array[Int] cluster_distances = [20, 10, 5] # CHANGING THIS WILL BREAK SECOND SCRIPT!
	String arg_denylist = if defined(persistent_denylist) then "--dl ~{persistent_denylist}" else ""
	String arg_shareemail = if defined(shareemail) then "-s ~{shareemail}" else ""
	String arg_microreact = if upload_clusters_to_microreact then "--upload_to_microreact" else ""
	String arg_disable_dropped_sample_failsafe = if no_dropped_sample_failsafe then "--no_dropped_sample_failsafe" else ""
	String arg_force_mr_update = if force_microreact_update then "--force_mr_update" else ""
	String arg_verbose = if verbose then "--verbose" else ""

	# naturally, this doesn't work on Cromwell
	#String? microreact_columns_csv = if defined(microreact_metadata_columns) then sep(",", microreact_metadata_columns) else ""
	
	command <<<
		set -eux pipefail
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting task"

		MICROREACT_COLUMNS_CSV=~{sep="," microreact_metadata_columns}

		if [[ "~{use_hardcoded_column_renames}" = "true" && "~{sample_metadata_tsv}" != "" ]]
		then
			if [[ $MICROREACT_COLUMNS_CSV = "Epi_Duplication,Year_Collected,Patient_County,State,Country,tbd_strain_per_tbprof,tbd_resistance,Submitter_Facility,Submitter_Facility_Sample_ID,Sequencing_Facility,Latitude,Longitude" ]]
			then
				sed -i '1s/tbd_strain_per_tbprof/Lineage_mycoTBProfiler/g; 1s/tbd_resistance/Resistance_mycoTBProfiler/g' "~{sample_metadata_tsv}"
			else
				echo "Your Microreact columns string: $MICROREACT_COLUMNS_CSV"
				echo "Expected: Epi_Duplication,Year_Collected,Patient_County,State,Country,tbd_strain_per_tbprof,tbd_resistance,Submitter_Facility,Submitter_Facility_Sample_ID,Sequencing_Facility,Latitude,Longitude"
				echo "You set use_hardcoded_column_renames to True, but your Microreact columns string doesn't match what we expect!"
				echo "Column renames are be hardcoded for safety due to a Cromwell bug (https://github.com/broadinstitute/cromwell/issues/7883). Crashing!"
				exit 1
			fi
		fi
			
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

			if [[ -f "~{microreact_update_template_json}" ]]
			then
				MR_UPDATE_JSON_ARG="--mr_update_template ~{microreact_update_template_json}"
			else
				echo "Upload to microreact is true, but no microreact_update_template_json provided. Crashing!"
				exit 1
			fi

			if [[ -f "~{microreact_blank_template_json}" ]]
			then
				MR_BLANK_JSON_ARG="--mr_blank_template ~{microreact_blank_template_json}"
			else
				echo "Upload to microreact is true, but no microreact_blank_template_json provided. Crashing!"
				exit 1
			fi

			if [[ -f "~{microreact_decimated_template_json}" ]]
			then
				MR_DECIMATED_JSON_ARG="--mr_decimated_template ~{microreact_decimated_template_json}"
			else
				echo -n "Upload to microreact is true, but no microreact_decimated_template_json provided. This isn't recommended,"
				echo -n "because decimated clusters on Microreact will never be updated, which might lead to incorrect assumptions."
			fi

			if [[ -f "~{cluster_subtrees_randomIDs_tarball}" && -f "~{cluster_matrices_randomIDs_tarball}" ]]
			then
				echo "Found cluster subtree and matrix tarballs; can upload to Microreact"
			else
				echo -n "Upload to microreact is true, but either cluster_subtrees_randomIDs_tarball or cluster_matrices_randomIDs_tarball "
				echo -n "(or both) is missing. Although technically optional if just finding persistent IDs, these files are necessary to show "
				echo -n "subtrees and distance matrices on Microreact. Crashing!"
				exit 1
			fi
		else
			TOKEN_ARG=""
			MR_UPDATE_JSON_ARG=""
			MR_BLANK_JSON_ARG=""
			MR_DECIMATED_JSON_ARG=""
		fi

		# we do similar logic within process_clusters.py too, but if we can crash before find_clusters.py that'd be ideal
		if [[ -f "~{persistent_ids}" ]]
		then
			if [[ -f "~{persistent_cluster_meta}" ]]
			then
				PERSISTENTIDS_ARG="--persistentids ~{persistent_ids}"
				PERSISTENTMETA_ARG="--persistentclustermeta ~{persistent_cluster_meta}"
			else
				echo "Found persistent IDs file but no persistent cluster meta. You need neither or both. Crashing!"
				exit 1
			fi
		else
			if [[ -f "~{persistent_cluster_meta}" ]]
			then
				echo "Found persistent cluster meta file but no persistent IDs. You need neither or both. Crashing!"
				exit 1
			else
				echo "Found neither persistent IDs file nor persistent cluster meta, will be running without persistent IDs"
				PERSISTENTIDS_ARG=""
				PERSISTENTMETA_ARG=""
			fi
		fi

		if [[ -f "~{sample_metadata_tsv}" ]]
		then
			SAMPLEMETADATA_ARG="--samplemeta ~{sample_metadata_tsv}"
		else
			SAMPLEMETADATA_ARG=""
		fi

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting newick (A_big.nwk) from input PB"
		matUtils extract -i ~{input_mat_with_new_samples} -t A_big.nwk
		cp ~{input_mat_with_new_samples} .

		if [[ -f "~{cluster_matrices_randomIDs_tarball}" ]]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expanding cluster matrix tarball..."
			pigz -dc "~{cluster_matrices_randomIDs_tarball}" | tar xf -
		fi
		
		if [[ -f "~{cluster_subtrees_randomIDs_tarball}" ]]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expanding cluster subtree tarball..."
			pigz -dc "~{cluster_subtrees_randomIDs_tarball}" | tar xf -
		fi

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

		if [[ ! "~{override_mass_rename_script}" == '' ]]
		then
			touch /HOME/ash/scripts/mass_rename_to_persistent_id.py
			rm /HOME/ash/scripts/mass_rename_to_persistent_id.py
			mv "~{override_mass_rename_script}" /HOME/ash/scripts/mass_rename_to_persistent_id.py
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

		mkdir logs

		echo "Contents of workdir:"
		tree

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generated these args for process_clusters.py:"
		echo "--combineddiff ~{combined_diff_file}"
		echo "--latestsamples ~{latest_samples_tsv}"
		echo "--latestclustermeta ~{latest_clusters_tsv}"
		echo "--mat_tree ~{input_mat_with_new_samples}"
		echo "--today ~{datestamp}"
		echo "~{arg_denylist}"
		echo "~{arg_disable_dropped_sample_failsafe}"
		echo "~{arg_verbose}"
		echo "$PERSISTENTIDS_ARG $PERSISTENTMETA_ARG"
		echo "$SAMPLEMETADATA_ARG"
		
		# microreact stuff
		echo "--mr_metadata_columns $MICROREACT_COLUMNS_CSV"
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
			--latestsamples "~{latest_samples_tsv}" \
			--latestclustermeta "~{latest_clusters_tsv}" \
			--mat_tree "~{input_mat_with_new_samples}" \
			--today ~{datestamp} \
			~{arg_denylist} \
			~{arg_disable_dropped_sample_failsafe} \
			~{arg_verbose} \
			$PERSISTENTIDS_ARG $PERSISTENTMETA_ARG \
			$SAMPLEMETADATA_ARG \
			--mr_metadata_columns $MICROREACT_COLUMNS_CSV \
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

		if [ ~{verbose} = "true" ]; then tree; fi

		# if process_clusters.py errored, NOW we should crash, since we have logs and such
		if [ "$PY_EXIT_CODE" -ne 0 ]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Crashing with rc $PY_EXIT_CODE because that's what process_clusters.py returned"
			exit $PY_EXIT_CODE
		fi

		if [ "~{previous_run_cluster_json}" != "" ]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running summarize_changes_alt.py"
			python3 /HOME/ash/scripts/summarize_changes_alt.py "all_cluster_information~{datestamp}.json"
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished summarize_changes_alt.py"
		fi

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running mass_rename_to_persistent_id.py"
		python3 /HOME/ash/scripts/mass_rename_to_persistent_id.py "~{arg_verbose}" --json "all_cluster_information~{datestamp}.json"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished mass_rename_to_persistent_id.py"

		echo "The IDs of these clusters were processed by process_clusters.py on ~{datestamp}(ish) and DO account for persistent cluster IDs. " > readme.txt
		echo "Note datestamp is set at execution of first task to ensure all outputs have same datestamp per workflow run, hence -ish. " >> readme.txt

		if [ ~{verbose} = "true" ]; then tree; fi

		find . -maxdepth 1 \( -name "a*.nwk" -o -name "a*.pb" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > "persisID_cluster_trees~{datestamp}.tar.gz"
		find . -maxdepth 1 \( -name "b*.nwk" -o -name "b*.pb" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > "persisID_cluster_trees_backmasked~{datestamp}.tar.gz"
		find . -maxdepth 1 \( -name "a*_dmtrx" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > "persisID_cluster_matrices~{datestamp}.tar.gz"
		find . -maxdepth 1 \( -name "b*_dmtrx" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > "persisID_cluster_matrices_backmasked~{datestamp}.tar.gz"

		# shellcheck disable=SC2317
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished task"

	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev14"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		# The amount of outputs we originally had was overloading Terra, so some of these are commented out now.
		# Also, we try to avoid globbing where possible to make finding outs in Terra bucket easier since globs
		# create a folder with a randomized name, which is annoying!

		###### IMPORTANT FILES THAT SHOULD ALWAYS GO INTO SUBSEQUENT RUNS IF THEY EXIST ######
		File new_persistent_ids = "persistentIDS" + datestamp + ".tsv"
		File new_persistent_meta = "persistentMETA" + datestamp + ".tsv"
		File final_cluster_information_json = "all_cluster_information" + datestamp + ".json"

		# cluster-specific subtrees and nwks
		# there is no internally masked big tree because masking is done per-cluster
		File cluster_trees_persisIDs = "persisID_cluster_trees" + datestamp + ".tar.gz"                        
		File cluster_trees_persisIDs_backmasked = "persisID_cluster_trees_backmasked" + datestamp + ".tar.gz"
		File cluster_matrices_persisIDs = "persisID_cluster_matrices" + datestamp + ".tar.gz"
		File cluster_matrices_persisIDs_backmasked = "persisID_cluster_matrices_backmasked" + datestamp + ".tar.gz"

		# stats for the nerds, currently unused downstream
		File new_samples_cluster_information = "new_samples" + datestamp + ".tsv"
		File all_samples_cluster_information = "all_samples" + datestamp + ".tsv"
		File? decimated_clusters = "decimated" + datestamp + ".tsv"
		File? updated_mr_URIs_file = "updated_mr_URIs" + datestamp + ".txt"

		# debug
		File? logs = "logs.zip"
		File? change_report_json       = "change_report" + datestamp + ".json"
		File? change_report_full       = "change_report_full"+datestamp+".txt"  # all clusters
		File? change_report_cdph       = "change_report_cdph"+datestamp+".txt"  # excludes 20-clusters
		#File? input_metadata_tsv       = sample_metadata_tsv                    # because metadata is mutable on Terra

		# can be used to annotate the Big Tree by cluster; uses Nextstrain (Auspice) metadata format
		# we currently don't annotate with Nextstrain as CDPH prefers Microreact, but the files are here if you want 'em
		File? samp_cluster_twn = "samp_persis20cluster" + datestamp + ".tsv"
		File? samp_cluster_ten = "samp_persis10cluster" + datestamp + ".tsv"
		File? samp_cluster_fiv = "samp_persis5cluster"  + datestamp + ".tsv"
	}
}