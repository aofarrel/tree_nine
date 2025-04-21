



# maybe put something like this after persistent cluster perl script to ensure that the IDs of merge/splits dont get reused



# ensure that we didn't have a total swap of cluster names due to samples being removed (kind of)
    # this is very tricky because we do not want to do what the persistent ID perl script is doing
    temp_latest_groupby_cluster = all_latest_samples.group_by("latest_cluster_id", maintain_order=True).agg([pl.col("sample_id"), pl.col("cluster_distance").unique()]).rename({"latest_cluster_id": "cluster_id"})
    persis_groupby_cluster = all_persistent_samples.group_by("cluster_id", maintain_order=True).agg(pl.col("sample_id"), pl.col("cluster_distance").unique().first())
    print_df_to_debug_log("Latest IDs, grouped by cluster", temp_latest_groupby_cluster)
    print_df_to_debug_log("Persistent IDs, grouped by cluster", persis_groupby_cluster)
    existing_cluster_ids = set(temp_latest_groupby_cluster["cluster_id"].to_list()) | set(persis_groupby_cluster["cluster_id"].to_list())
    latest_overrides = {}
    all_latest_samples_set = set(all_latest_samples["latest_cluster_id"].to_list())

for latest_row in temp_latest_groupby_cluster.iter_rows(named=True):
        latest_cluster_id = latest_row["cluster_id"]
        latest_cluster_distance = latest_row["cluster_distance"][0]
        this_clusters_latest_samps = set(latest_row["sample_id"])
        this_clusters_persis_row = persis_groupby_cluster.filter(pl.col("cluster_id") == latest_cluster_id)

        if this_clusters_persis_row.is_empty():
            logging.warning("%s present in latest but not persistent", latest_cluster_id)
        else:
            this_clusters_persis_samps = set(this_clusters_persis_row["sample_id"].item())
            # In and of itself, a disjoint union could be okay. But we have to watch out for decimated
            # clusters that lose all (or all except one, since a "cluster" of 1 is no longer a cluster)
            # of their samples, because we don't want to risk reusing their IDs.
            # Keep in mind that "all" latest samples is all latest samples that cluster. It excludes
            # unclustered samples.
            # Previously we checked if
            # this_clusters_latest_samps.isdisjoint(persis_groupby_cluster.filter(pl.col("cluster_id") == latest_cluster_id)),
            # but this was a terrible idea because it would fire whenever samples simply generated a different cluster ID.
            # We instead want to ask "are there samples "
            if this_clusters_latest_samps.isdisjoint(this_clusters_persis_samps):
                logging.debug("Disjoint union of %s", latest_cluster_id)
                #logging.debug("    This run, %s had these samples: %s", latest_cluster_id, this_clusters_latest_samps)
                #logging.debug("    But previous run, %s had these samples: %s", latest_cluster_id, this_clusters_persis_samps)
                intersection = all_latest_samples_set.intersection(this_clusters_persis_samps)
                if len(intersection) == 0:
                    logging.warning("Persistent cluster %s@%s with samples %s is decimated", 
                        latest_cluster_id, latest_cluster_distance, this_clusters_persis_samps)
                    new_cluster_id = generate_truly_unique_cluster_id(existing_cluster_ids, args.denylist)
                    latest_overrides[cluster_id] = new_cluster_id
                    existing_cluster_ids.add(new_cluster_id)

                    # TODO: change back to check = True once done debugging
                    logging.warning("Generated new workdir cluster ID: %s → %s", latest_cluster_id, new_cluster_id)
                    subprocess.run(f"mv a{latest_cluster_id}.nwk a{new_cluster_id}.nwk", shell=True, check=False)
                    logging.info("Renamed nwk")
                    subprocess.run(f"mv a{latest_cluster_id}.pb a{new_cluster_id}.pb", shell=True, check=False)
                    logging.info("Renamed pb")
                    subprocess.run(f"mv a{latest_cluster_id}_dmtrx.tsv a{new_cluster_id}_dmtrx.tsv", shell=True, check=False)
                    logging.info("Renamed distance matrix")
                elif len(intersection) == 1:
                    logging.warning("Persistent cluster %s@%s with samples %s is near-decimated; sole remaining sample is %s", 
                        cluster_id, cluster_distance, this_clusters_persis_samps, golden_sample)
                    # WARNING: UNTESTED!!!
                    # If this is a cluster that loses all of its samples except one, and no new samples cluster with the
                    # one remaining sample, we consider it decimated. If at least one new sample clusters to it at the 
                    # same distance, it is not decimated. But as noted above, "all" latest samples is all latest samples
                    # that cluster, so if intersection is len(1), then that 1 sample has clustered...somewhere. The issue is
                    # that it might be a cluster of a different distance.
                    # For example: A, B, C formed a persistent 20 SNP cluster 000001, while A and B formed persistent 10 SNP
                    # cluster 000002. In the latest run, B and C are removed, and A clusters with new sample D in a 20 SNP
                    # cluster currently called 000050. all_latest_samples_set.intersection(this_clusters_persis_samps) = set(A)
                    # for 000050, so Marc's script should translate 000050 -> 000001. But when we inner join the persistent
                    # 10-clusters with the latest 10-clusters, A is not in a 10 cluster, so it's excluded and not factored in.
                    # This could lead to improper reusage of cluster ID 000002.
                    golden_sample = intersection.pop()
                    latest_samples_of_this_distance = all_latest_samples.filter(pl.col("cluster_distance") == cluster_distance)
                    if golden_sample in latest_samples_of_this_distance["cluster_id"].to_list():
                        # do nothing except throw warning
                        logging.warning("However, it has picked up new samples at the correct distance.")
                    else:
                        logging.warning("It lacks samples at the correct distance and will be treated as if it were decimated.")

                        # TODO: change back to check = True once done debugging
                        logging.warning("Generated new workdir cluster ID: %s → %s", latest_cluster_id, new_cluster_id)
                        subprocess.run(f"mv a{latest_cluster_id}.nwk a{new_cluster_id}.nwk", shell=True, check=False)
                        logging.info("Renamed nwk")
                        subprocess.run(f"mv a{latest_cluster_id}.pb a{new_cluster_id}.pb", shell=True, check=False)
                        logging.info("Renamed pb")
                        subprocess.run(f"mv a{latest_cluster_id}_dmtrx.tsv a{new_cluster_id}_dmtrx.tsv", shell=True, check=False)
                        logging.info("Renamed distance matrix")
    
    if len(latest_overrides) > 0:
        all_latest_samples = all_latest_samples.with_columns(
            pl.when(pl.col("latest_cluster_id").is_in(list(latest_overrides.keys())))
            .then(pl.col("latest_cluster_id").replace(latest_overrides))
            .otherwise(pl.col("latest_cluster_id"))
            .alias("latest_cluster_id")
        )
        print_df_to_debug_log("Latest IDs after name changes", all_latest_samples)

        # write to denylist so subsequent runs will still avoid this cluster name
        if args.denylist is not None:
            with open(args.denylist, "a") as denylist:
                for key in latest_overrides:
                    denylist.write(str(key) + '\n')
        else:
            with open('clusterid_denylist.txt', 'w') as denylist:
                for key in latest_overrides:
                    denylist.write(str(key) + '\n')
    else:
        logging.debug("Did not change any persistent ID names prior to running the main script that also changes persistent IDs (just roll with it)")