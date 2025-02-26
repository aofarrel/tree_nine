# FIND CLUSTERS - VERSION 1.6.1
SCRIPT_PATH = '/scripts/find_clusters.py'


# pylint: disable=too-complex,too-many-locals,too-many-statements
# pylint: disable=multiple-statements,consider-using-sys-exit,wrong-import-position,no-else-return,unnecessary-pass,consider-using-enumerate,useless-suppression

import os
import argparse
import logging
from datetime import date
import time
import subprocess
import numpy as np
import ete3
import pandas as pd # not messing with polars' restrictions on TSVs today no sir

INT16_MAX = np.iinfo(np.int16).max
INT32_MAX = np.iinfo(np.int32).max
INT64_MAX = np.iinfo(np.int64).max

def main():
    parser = argparse.ArgumentParser(description="Clusterf...inder")
    parser.add_argument('mat_tree', type=str, required=True, help='input MAT (.pb)')
    parser.add_argument('nwk_tree', type=str, required=True, help='input nwk (.nwk)')
    parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')
    parser.add_argument('-d', '--distance', default=20, type=int, help='max distance between samples to identify as clustered')
    parser.add_argument('-rd', '--recursive-distance', type=lambda x: [int(i) for i in x.strip('"').split(',')], help='after identifying --distance cluster, search for subclusters with these distances')
    parser.add_argument('-c', '--contextsamples', default=0, type=int, help='number of context samples to add per cluster (appears only in nwk)')
    parser.add_argument('-t', '--type', choices=['BM', 'NB'], type=str.upper, help='BM=backmasked, NB=not-backmasked; will add BM/NB before prefix')
    parser.add_argument('-cn', '--collection-name', default='unnamed', type=str, help='name of this group of samples (do not include a/b prefix)')
    parser.add_argument('-sf', '--startfrom', default=0, type=int, help='the six-digit int part of cluster UUIDs will begin with the next integer after this one')
    parser.add_argument('-v', '--verbose', action='store_true', help='enable info logging')
    parser.add_argument('-vv', '--veryverbose', action='store_true', help='enable debug logging')

    # mostly used for recursion
    parser.add_argument('-nb', '--nobackmask', action='store_true', help='do not backmask')
    parser.add_argument('-nm', '--nomatrix', action='store_true', help='do not matrix (if clustering, clusters will still have their own matrices written)')
    parser.add_argument('-nc', '--nocluster', action='store_true', help='do not search for clusters (will not recurse)')
    parser.add_argument('-nl', '--nolonely', action='store_true', help='do not make a subtree for unclustered samples')
    parser.add_argument('-neo', '--noextraouts', action='store_true', help='do not write extra summary information to their own file')

    args = parser.parse_args()
    type_prefix, args_dot_type, matrix_out, tree_out, clus_distance_i32 = initialize(args)

    t = ete3.Tree(args.nwk_tree, format=1)
    logging.debug("Input nwk_tree loaded")
    samps = args.samples.split(',') if args.samples else sorted([leaf.name for leaf in t])
    samps, matrix, clusters, lonely = dist_matrix(t, samps, args, clus_distance_i32)
    if args.nocluster:
        logging.debug("Processed %s samples (not clustering)", len(samps)) # sort of a misnomer since we already clustered...
    else:
        logging.info("Processed %s samples, found %s clusters @ %s SNPs", len(samps), len(clusters), clus_distance_i32)
    #logging.debug("Samples processed: %s", samps) # check if alphabetized

    # write distance matrix (and subtree)
    # this can be called during recursion -- it may be a dmatrix/tree of the overall sample pool, or a cluster/subcluster
    if not args.nomatrix:
        write_matrix(samps, matrix, matrix_out)
        write_subtree(samps, args.mat_tree, tree_out, args.collection_name)        

    # backmask (and write BM matrix) if possible
    if args.type != 'BM' and not args.nobackmask:
#            bm_samps = ["[BM]" + s for s in samps]
#            samples_in_cluster_str = ','.join(bm_samps)
            # if we are in a (sub)cluster recursion, we know that the UUID of this (sub)cluster will be args.startfrom minus 1
            # if we are backmasking the entire tree, we assume args.startfrom is not defined. is this iffy? mmmm yeah but it'll do for now
            if args.startfrom: # "b" will be added by the called script
                this_backmasked_groups_name = str(int(args.startfrom) - 1)
            else:
                this_backmasked_groups_name = "tree"
            handle_subprocess(f"Backmasking {this_backmasked_groups_name} group of samples",
                f"python3 ./scripts/diffdiff.py -b --mask_outfile {matrix_out}")
#            # We now have backmasked diff files with extension .backmask.diff, and within those files, their sample names
#            # have changed to include [BM]
#
#            logging.debug("Backmask group name: %s", this_backmasked_groups_name)
#            logging.debug("Backmasked samples in cluster str: %s", samples_in_cluster_str)
#
#            handle_subprocess("Generating a backmasked distance matrix and subcluster tree...", # TODO: do I need a wholeahh new mat tree too?
#                f"python3 {SCRIPT_PATH} '{args.mat_tree}' '{this_backmasked_groups_name}' -d {clus_distance_i32} -s{samples_in_cluster_str} -v -t BM --collection-name {this_backmasked_groups_name} -neo -nl -nb")
#            samples_in_cluster_str = None


    # this could probably be made more efficient, but it's not worth refactoring
    if not args.nocluster:
        # sample_cluster is the Nextstrain-style TSV used for annotation, eg:
        # sample12    cluster1
        # sample13    cluster1
        # sample14    cluster1
        sample_cluster = ['Sample\tCluster\n']
        sample_clusterUUID = ['Sample\tClusterUUID\n']

        # cluster_samples is for matutils extract to generate subtrees, eg:
        # cluster1    sample12,sample13,sample14
        cluster_samples = ['Cluster\tSamples\n']

        # summary information for humans to look at
        n_clusters = len(clusters) # immutable
        n_samples_in_clusters = 0  # mutable

        number_start = args.startfrom # mutable if there's subclusters
        number_start_old = number_start

        for n in range(n_clusters):
            # get basic information -- we can safely sort here as do not use the array directly
            samples_in_cluster = sorted(list(clusters[n]))
            assert len(samples_in_cluster) == len(set(samples_in_cluster))
            n_samples_in_clusters += len(samples_in_cluster) # samples in ANY cluster, not just this one
            samples_in_cluster_str = ",".join(samples_in_cluster)

            # TODO: this heuristic will eventually be replaced by metadata
            is_cdph = is_cdph = any(
                samp_name[:2].isdigit() or
                (samp_name.startswith("[BM]") and samp_name[4:6].isdigit()) or
                (samp_name.startswith("[NB]") and samp_name[4:6].isdigit())
                for samp_name in samples_in_cluster
            )
            locale = 'CA' if is_cdph else '??'
            #number_part = n + number_start
            number_part = number_start
            UUID = f"{str(clus_distance_i32).zfill(2)}SNP-{locale}-{str(date.today().year)}M{str(date.today().month).zfill(2)}-{str(number_part).zfill(6)}"
            cluster_name = UUID # previously they were meaningfully different, but this is less prone to nonsense
            logging.info("Identified %s with %s members", cluster_name, len(samples_in_cluster))

            # build cluster_samples line for this cluster
            cluster_samples.append(f"{cluster_name}\t{samples_in_cluster_str}\n")

            # build sample_cluster lines for this cluster - this will be used for auspice (JSON) annotation
            # TODO: this gets cleared per recursion!!!!
            for s in samples_in_cluster:
                sample_cluster.append(f"{s}\t{cluster_name}\n")
                sample_clusterUUID.append(f"{s}\t{UUID}\n")

            # this is an optional thing for matUtils extract, calculate now might as well sure
            minimum_tree_size = args.contextsamples + len(samples_in_cluster)

            # write as much information about this cluster as we have
            current_cluster_headers = 'latest_cluster_id\tcurrent_date\tcluster_distance\tn_samples\tminimum_tree_size\tsample_ids\n'
            if not os.path.isfile("latest_clusters.tsv"):
                with open("latest_clusters.tsv", "w", encoding="utf-8") as current_clusters:
                    current_clusters.write(current_cluster_headers)
            with open("latest_clusters.tsv", "a", encoding="utf-8") as current_clusters:  # TODO: eventually add old/new samp information
                current_clusters.write(f"{cluster_name}\t{date.today().isoformat()}\t{clus_distance_i32}\t{len(samples_in_cluster)}\t{minimum_tree_size}\t{samples_in_cluster_str}\n")

            # do something similar for samples
            current_sample_headers = 'sample_id\tcluster_distance\tlatest_cluster_id\n'
            if not os.path.isfile("latest_samples.tsv"):
                with open("latest_samples.tsv", "w", encoding="utf-8") as latest_samples:
                    latest_samples.write(current_sample_headers)
            with open("latest_samples.tsv", "a", encoding="utf-8") as latest_samples:  # TODO: eventually add old/new samp information
                for s in samples_in_cluster:
                    latest_samples.write(f"{s}\t{clus_distance_i32}\t{cluster_name}\n")

            # generate the distance matrix for the CURRENT cluster
            handle_subprocess(f"Generating {cluster_name}'s distance matrix...",  f"python3 {SCRIPT_PATH} '{args.mat_tree}' '{args.nwk_tree}' -d {clus_distance_i32}  -s{samples_in_cluster_str} -v {args_dot_type} -neo -nc -nl -cn {cluster_name}")
            time.sleep(1)

            # recurse as needed for subclusters
            logging.debug('args.recursive_distance: %s', args.recursive_distance)
            if args.recursive_distance is None or len(args.recursive_distance) == 0:
                number_start = number_start + 1
            elif len(args.recursive_distance) == 1:
                next_recursion = args.recursive_distance[0]
                logging.info("Looking for %s-SNP subclusters...", next_recursion)
                result = subprocess.run(f"python3 {SCRIPT_PATH} '{args.mat_tree}' '{args.nwk_tree}' -d {next_recursion} -s{samples_in_cluster_str} -v {args_dot_type} -neo -nm -nl -sf {number_part+1}", shell=True, check=False)
                number_start = result.returncode
            elif len(args.recursive_distance) == 2:
                next_recursion = args.recursive_distance[0]
                next_next_recursion = args.recursive_distance[1]
                logging.info("Looking for %s-SNP subclusters...", next_recursion)
                result = subprocess.run(f"python3 {SCRIPT_PATH} '{args.mat_tree}' '{args.nwk_tree}' -d {next_recursion} -s{samples_in_cluster_str} -v {args_dot_type} -neo -nm -nl -sf {number_part+1} -rd {next_next_recursion}", shell=True, check=False)
                number_start = result.returncode
            else:
                subsequent_recursions = f'-rd {",".join(map(str, args.recursive_distance[:1]))}' # pylint: disable=bad-builtin
                logging.info("Looking for %s-SNP subclusters...", args.recursive_distance[0])
                result = subprocess.run(f"python3 {SCRIPT_PATH} '{args.mat_tree}' '{args.nwk_tree}' -d {next_recursion} -s{samples_in_cluster_str} -v {args_dot_type} -neo -nm -nl -sf {number_part+1} {subsequent_recursions}", shell=True, check=False)
                number_start = result.returncode
            logging.debug('%s done recursing, number_start was %s, number_part is %s, number_start now %s', cluster_name, number_start_old, number_part, number_start)

        # add in the unclustered samples (outside for loop to avoid writing multiple times)
        # however, don't add to the UUID list, or else persistent cluster IDs will break
        if not args.nolonely:
            lonely = sorted(list(lonely))
            for george in lonely: # W0621, https://en.wikipedia.org/wiki/Lonesome_George
                sample_cluster.append(f"{george}\tlonely\n")
            with open(f"{type_prefix}_lonely.txt", "w", encoding="utf-8") as unclustered_samples_list:
                unclustered_samples_list.writelines(line + '\n' for line in lonely)
            cluster_samples.append(f"lonely\t{','.join(lonely)}\n")
            lonely_minimum_tree_size = int(args.contextsamples) + len(lonely)
            handle_subprocess("Also extracting a tree for lonely samples...", # we'll accept -N here as subsubtrees might actually be helpful here
                f'matUtils extract -i "{args.mat_tree}" -t "LONELY" -s {type_prefix}_lonely.txt -N {lonely_minimum_tree_size}')
            os.rename("subtree-assignments.tsv", "lonely-subtree-assignments.tsv")
            [os.rename(f, f[:-2] + "nwk") for f in os.listdir() if f.endswith(".nw")] # pylint: disable=expression-not-assigned
            logging.debug(os.listdir('.'))

        if not args.noextraouts:
            logging.info("Writing auspice-style TSV for annotation of clusters...")
            with open("temp_cluster_anno.tsv", "a", encoding="utf-8") as samples_for_annotation:
                samples_for_annotation.writelines(sample_cluster)

            # Because we are recursing, samples can have more than one subtree assignment. As such, we don't use cluster_samples
            # for usher extraction anymore.
            # We also don't bother with sample_clusterUUID either as we don't need a file that excludes unclustered samples;
            # we already have latest_samples.tsv for that!

            logging.info("Writing other little doodads...")
            with open("n_big_clusters", "w", encoding="utf-8") as n_cluster: n_cluster.write(str(n_clusters))
            with open("n_samples_in_clusters", "w", encoding="utf-8") as n_cluded: n_cluded.write(str(n_samples_in_clusters))
            with open("n_samples_processed", "w", encoding="utf-8") as n_processed: n_processed.write(str(len(samps)))
            with open("n_unclustered", "w", encoding="utf-8") as n_lonely: n_lonely.write(str(len(lonely)))
            find_neighbors(matrix, samps, f"{args.collection_name}_neighbors.tsv")
    else:
        logging.debug("Not clustering, exiting...")

    if not args.nocluster:
        last_cluster_number = number_start
        logging.debug('returning %s @ %sSNP', last_cluster_number, clus_distance_i32)
    else:
        last_cluster_number = 99999999999999999999999999999999999999999999999999999999999  # should never get read, but if it does, we'll certainly notice!
    exit(last_cluster_number)

def initialize(args):
    logging.basicConfig(level=logging.DEBUG if args.veryverbose else logging.INFO if args.verbose else logging.WARNING)
    if args.type == 'BM':
        type_prefix, args_dot_type = 'b', '-t BM' # for "backmasked"
    elif args.type == 'NB':
        type_prefix, args_dot_type = 'a', '-t NB' # for... uh... Absolutelynotbackmasked
    else:
        type_prefix, args_dot_type = '', ''
    matrix_out = f"{type_prefix}{args.collection_name}_dmtrx.tsv"
    tree_out = f"{type_prefix}{args.collection_name}" # note that extension breaks if using -N, see https://github.com/yatisht/usher/issues/389
    if args.distance > INT32_MAX:
        logging.error("args.distance is a value greater than the signed-int32 maximum used when generating matrices; cannot continue")
        exit(1)
    else:
        clus_distance_i32 = np.int32(args.distance)
    logging.info(f"~~ Processing {type_prefix}{args.collection_name}#N+{args.startfrom} @ {clus_distance_i32} ~~") # pylint: disable=logging-fstring-interpolation
    return type_prefix, args_dot_type, matrix_out, tree_out, clus_distance_i32

def write_matrix(samps, matrix, matrix_out):
    # Write distance matrix
    assert not os.path.exists(matrix_out), f"Tried to write {matrix_out} but it already exists?!"
    with open(matrix_out, "a", encoding="utf-8") as outfile:
        outfile.write('sample\t'+'\t'.join(samps))
        outfile.write("\n")
        for k in range(len(samps)): # don't change to enumerate without changing i; with enumerate it's a tuple
            #strng = np.array2string(matrix[i], separator='\t')[1:-1]
            line = [ str(int(count)) for count in matrix[k]]
            outfile.write(f'{samps[k]}\t' + '\t'.join(line) + '\n')
        time.sleep(1)
    logging.debug("Wrote distance matrix: %s", matrix_out)

def write_subtree(samps, input_mat_tree, tree_out, collection_name):
    # TODO: also extract JSON version of the tree and add metadata to it (-M metadata_tsv) even though that doesn't go to MR
    #minimum_tree_size = args.contextsamples + len(samps) # this variable is overwritten in cluster loop # THIS CAUSES SUB-SUB-TREES
    assert not os.path.exists(f"{tree_out}.nw"), f"Tried to make subtree called {tree_out}.nw but it already exists?!"
    with open("temp_extract_these_samps.txt", "w", encoding="utf-8") as temp_extract_these_samps:
        temp_extract_these_samps.writelines(line + '\n' for line in samps)
     handle_subprocess(f"Extracting {tree_out} pb for {collection_name}...",
        f'matUtils extract -i "{input_mat_tree}" -o "{tree_out}.pb" -s temp_extract_these_samps.txt') # removed -N {minimum_tree_size}
    handle_subprocess(f"Extracting {tree_out} nwk for {collection_name}...",
        f'matUtils extract -i "{input_mat_tree}" -t "{tree_out}.nwk" -s temp_extract_these_samps.txt') # removed -N {minimum_tree_size}
    if os.path.exists(f"{tree_out}-subtree-1.nw"):
        logging.warning("Generated multiple subtrees for %s, attempting batch rename (this may break things)", collection_name)
        [os.rename(f, f[:-2] + "nwk") for f in os.listdir() if f.endswith(".nw")] # pylint: disable=expression-not-assigned
    else:
        [os.rename(f, f[:-13] + ".nwk") for f in os.listdir() if f.endswith("-subtree-0.nw")] # pylint: disable=expression-not-assigned
    if os.path.exists("subtree-assignments.tsv"):
        os.rename("subtree-assignments.tsv", "lonely-subtree-assignments.tsv")

def handle_subprocess(explainer, system_call_as_string):
    # Wrapper function for os.system() calls which uses time.sleep() to avoid some issues with parallelism
    logging.info(explainer)
    logging.debug(system_call_as_string)
    os.system(system_call_as_string)
    time.sleep(1)

def get_true_clusters(neighbors, args_no_cluster):
    # From neighbors we generated while making distance matrix, define clusters
    if not args_no_cluster:
        true_clusters, first_iter = [], True
        for pairs in neighbors:
            existing_cluster = False
            if first_iter:
                true_clusters.append(set([pairs[0], pairs[1]]))
            else:
                for sublist in true_clusters:
                    if pairs[0] in sublist:
                        sublist.add(pairs[1])
                        existing_cluster = True
                    if pairs[1] in sublist: # NOT ELSE IF
                        sublist.add(pairs[0])
                        existing_cluster = True
                if not existing_cluster:
                    true_clusters.append(set([pairs[0], pairs[1]]))
            first_iter = False
        return true_clusters
    else:
        return None

def path_to_root(ete_tree, node_name):
    # Browse the tree from a specific leaf to the root
    node = ete_tree.search_nodes(name=node_name)[0]
    path = [node]
    while node:
        node = node.up
        path.append(node)
    #logging.debug("path for %s: %s", node_name, path)
    return path

def overflow_check(this_path, that_path, this_samp, that_samp):
    # Checks for integer overflow
    total_distance = this_path + that_path
    if total_distance > INT32_MAX:
        logging.warning("Total distance between %s and %s is %s, greater than signed 32-bit maximum; will store as %s", this_samp, that_samp, total_distance, INT32_MAX)
        return np.int32(INT32_MAX)
    else:
        return np.int32(total_distance)


def dist_matrix(tree_to_matrix, samples, args, clus_distance_i32):
    # Generate a distance matrix as an np array, as well as a list of samples associated with it,
    # and which samples go into which cluster (if any)
    samp_ancs, neighbors, unclustered = {}, [], set()

    for sample in samples:
    #for sample in progressbar.tqdm(samples, desc="Finding roots and branch lengths"):
        s_ancs = path_to_root(tree_to_matrix, sample)
        samp_ancs[sample] = s_ancs

    # Currently using a 32-bit signed int matrix in hopes of less aggressive RAM usage
    matrix = np.full((len(samples),len(samples)), -1, dtype=np.int32)

    for i in range(len(samples)):
    #for i in progressbar.trange(len(samples), desc="Creating matrix"): # trange is a tqdm optimized version of range
        this_samp = samples[i]
        definitely_in_a_cluster = False

        for j in range(len(samples)):
            that_samp = samples[j]
            #Future goal: add catch to prevent reiteration of already checked pairs --> if matrix[i][j] == -1 should handle this
            if matrix[i][j] == -1: # ie, we haven't calculated this one yet
                if that_samp == this_samp: # self-to-self
                    matrix[i][j] = np.int32(0)
                else:
                    #find lca, add up branch lengths
                    this_path = 0
                    that_path = 0

                    for a in samp_ancs[this_samp]:
                        this_path += a.dist
                        if a in samp_ancs[that_samp]:
                            lca = a
                            this_path -= a.dist
                            break

                    for a in samp_ancs[that_samp]:
                        that_path += a.dist
                        if a == lca:
                            that_path -= a.dist
                            break

                    logging.debug("  sample %s vs other sample %s: this_path %s, that_path %s", this_samp, that_samp, this_path, that_path)
                    total_distance = overflow_check(this_path, that_path, this_samp, that_samp)
                    matrix[i][j], matrix[j][i] = total_distance, total_distance
                    if not args.nocluster and total_distance <= clus_distance_i32:
                        logging.debug("  %s and %s seem to be in a cluster (%s)", this_samp, that_samp, total_distance)
                        neighbors.append(tuple((this_samp, that_samp)))
                        definitely_in_a_cluster = True

        # after iterating through all of j, if this_sample (i) is not in a cluster, make note of that
        if not args.nocluster and not definitely_in_a_cluster:
            second_smallest_distance = np.partition(matrix[i], 1)[1] # second smallest, because smallest is self-self at 0
            if second_smallest_distance <= clus_distance_i32:
                logging.debug("  Oops, %s was already clustered! (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                pass
            else:
                logging.debug("  %s appears to be truly unclustered (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                unclustered.add(this_samp)

    # finished iterating, let's see what our clusters look like
    true_clusters = get_true_clusters(neighbors, args.nocluster) # None if args.noclusters
    logging.debug("Returning:\n\tsamples:\n%s\n\tmatrix:\n%s\n\ttrue_clusters:\n%s\n\tunclustered:\n%s", samples, matrix, true_clusters, unclustered)
    return samples, matrix, true_clusters, unclustered

def find_neighbors(distance_matrix: np.ndarray, sample_names: list, output_tsv: str):
    logging.info("Searching for closest and furthest neighbor samples...")
    rows = []
    for i, sample in enumerate(sample_names):
        distances = distance_matrix[i]
        distances[i] = 9999999 # exclude self-self from closest by temporarily setting to something goofy
        closest_dist = np.min(distances)
        distances[i] = 0 # fix self-self
        furthest_dist = np.max(distances)
        closest_neighbors = [sample_names[j] for j in np.where(distances == closest_dist)[0]]
        farthest_neighbors = [sample_names[j] for j in np.where(distances == furthest_dist)[0]]
        rows.append([sample, ", ".join(closest_neighbors), closest_dist, ", ".join(farthest_neighbors), furthest_dist])
    df = pd.DataFrame(rows, columns=["sample", "closest_neighbor(s)", "closest_distance", "furthest_sample(s)", "furthest_distance"])
    print(df)
    df.to_csv(output_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
