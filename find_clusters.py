print("VERSION 1.5.6")
script_path = '/scripts/find_clusters.py'

import os
import argparse
import logging
from datetime import date
import time
import numpy as np
import ete3
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Clusterf...inder")
    parser.add_argument('mat_tree', type=str, help='input MAT')
    parser.add_argument('nwk_tree', type=str, help='input nwk')
    parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')
    parser.add_argument('-d', '--distance', default=20, type=int, help='max distance between samples to identify as clustered')
    parser.add_argument('-rd', '--recursive-distance', type=lambda x: [int(i) for i in x.strip('"').split(',')], help='after identifying --distance cluster, search for subclusters with these distances')
    parser.add_argument('-c', '--contextsamples', default=0, type=int, help='number of context samples to add per cluster (appears only in nwk)')
    parser.add_argument('-mt', '--microreacttokenfile', type=str, help='!!!PATH!! to microreact token file')
    parser.add_argument('-t', '--type', choices=['BM', 'NB'], type=str.upper, help='BM=backmasked, NB=not-backmasked; will add BM/NB before prefix')
    parser.add_argument('-mo', '--matrixout', type=str, help='outname for the big matrix (sans ext)')
    parser.add_argument('-sf', '--startfrom', default=0, type=int, help='the six-digit int part of cluster UUIDs will begin with the next integer after this one')
    parser.add_argument('-v', '--verbose', action='store_true', help='enable info logging')
    parser.add_argument('-vv', '--veryverbose', action='store_true', help='enable debug logging')
    
    # mostly used for recursion
    parser.add_argument('-nm', '--nomatrix', action='store_true', help='do not matrix (if clustering, clusters will still have their own matrices written)')
    parser.add_argument('-nc', '--nocluster', action='store_true', help='do not search for clusters (will not recurse)')
    parser.add_argument('-nl', '--nolonely', action='store_true', help='do not make a subtree for unclustered samples')
    parser.add_argument('-neo', '--noextraouts', action='store_true', help='do not write extra summary information to their own file')

    args = parser.parse_args()
    if args.veryverbose:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)
    
    #logging.debug("~~Passed in~~")
    #logging.debug(args)
    #time.sleep(1)

    #logging.debug("Loading tree...")
    t = ete3.Tree(args.nwk_tree, format=1)
    #logging.debug("Tree loaded")
    samps = args.samples.split(',') if args.samples else sorted([leaf.name for leaf in t])
    if args.type == 'BM':
        is_backmasked, type_prefix, args_dot_type = True, 'BM_', '-t BM'
    elif args.type == 'NB':
        is_backmasked, type_prefix, args_dot_type = False, 'NB_', '-t NB'
    else:
        is_backmasked, type_prefix, args_dot_type = None, '', ''
    prefix = type_prefix
    matrix_out = args.matrixout if args.matrixout else f'{prefix}_unknown_matrix'

    samps, mat, clusters, lonely = dist_matrix(t, samps, args)
    total_samples_processed = len(samps)
    if args.nocluster:
        logging.debug("Processed %s samples (not clustering)", total_samples_processed)
    else:
        logging.info("Processed %s samples, found %s clusters @ %s SNPs", total_samples_processed, len(clusters), args.distance)
    #logging.debug("Samples processed: %s", samps) # check if alphabetized
    time.sleep(1)

    # write distance matrix
    if not args.nomatrix:
        with open(f"{matrix_out}.tsv", "a", encoding="utf-8") as outfile: # TODO: should this be in write mode instead of append?
            outfile.write('sample\t'+'\t'.join(samps))
            outfile.write("\n")
            for k in range(len(samps)): # don't change to enumerate without changing i; with enumerate it's a tuple
                #strng = np.array2string(mat[i], separator='\t')[1:-1]
                line = [ str(int(count)) for count in mat[k]]
                outfile.write(f'{samps[k]}\t' + '\t'.join(line) + '\n')
            time.sleep(1)
        logging.debug(f"Wrote distance matrix: {matrix_out}.tsv")

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
            is_cdph = is_cdph = any(
                samp_name[:2].isdigit() or 
                (samp_name.startswith("[BM]") and samp_name[4:6].isdigit()) or
                (samp_name.startswith("[NB]") and samp_name[4:6].isdigit())
                for samp_name in samples_in_cluster
            )
            locale = 'CA' if is_cdph else '??'
            if prefix == 'BM':
                short_prefix = 'm'
            elif prefix == 'NB':
                short_prefix = 'n'
            else:
                short_prefix = prefix
            #number_part = n + number_start
            number_part = number_start
            UUID = f"{str(args.distance).zfill(2)}SNP-CA-{str(date.today().year)}-{short_prefix}{str(number_part).zfill(6)}"
            cluster_name = UUID # previously they were meaningfully different, but this is less prone to nonsense
            logging.info("Identified %s with %s members", cluster_name, len(samples_in_cluster))

            # build cluster_samples line for this cluster
            cluster_samples.append(f"{cluster_name}\t{samples_in_cluster_str}\n")
            time.sleep(1)
            
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
                current_clusters.write(f"{cluster_name}\t{date.today().isoformat()}\t{args.distance}\t{len(samples_in_cluster)}\t{minimum_tree_size}\t{samples_in_cluster_str}\n")

            # do something similar for samples
            current_sample_headers = 'sample_id\tcluster_distance\tlatest_cluster_id\n'
            if not os.path.isfile("latest_samples.tsv"):
                with open("latest_samples.tsv", "w", encoding="utf-8") as latest_samples:
                    latest_samples.write(current_sample_headers)
            with open("latest_samples.tsv", "a", encoding="utf-8") as latest_samples:  # TODO: eventually add old/new samp information
                for s in samples_in_cluster:
                    latest_samples.write(f"{s}\t{args.distance}\t{cluster_name}\n")

            # generate the distance matrix for the CURRENT cluster
            handle_subprocess(f"Generating {cluster_name}'s distance matrix...",  f"python3 {script_path} '{args.mat_tree}' '{args.nwk_tree}' -d {args.distance}  -s{samples_in_cluster_str} -v {args_dot_type} -neo -nc -nl -mo {prefix}{cluster_name}_dmtrx")
            time.sleep(1)

            # recurse as needed for subclusters
            logging.debug('args.recursive_distance: %s', args.recursive_distance)
            if args.recursive_distance is None or len(args.recursive_distance) == 0:
                number_start = number_start + 1
            elif len(args.recursive_distance) == 1:
                next_recursion = args.recursive_distance[0]
                logging.info("Looking for %s-SNP subclusters...", next_recursion)
                result = subprocess.run(f"python3 {script_path} '{args.mat_tree}' '{args.nwk_tree}' -d {next_recursion} -s{samples_in_cluster_str} -v {args_dot_type} -neo -nm -nl -sf {number_part+1}", shell=True)
                number_start = result.returncode
            elif len(args.recursive_distance) == 2:
                next_recursion = args.recursive_distance[0]
                next_next_recursion = args.recursive_distance[1]
                logging.info(f"Looking for {next_recursion}-SNP subclusters...")
                result = subprocess.run(f"python3 {script_path} '{args.mat_tree}' '{args.nwk_tree}' -d {next_recursion} -s{samples_in_cluster_str} -v {args_dot_type} -neo -nm -nl -sf {number_part+1} -rd {next_next_recursion}", shell=True)
                number_start = result.returncode
            else:
                subsequent_recursions = f'-rd {",".join(map(str, args.recursive_distance[:1]))}'
                logging.info(f"Looking for {args.recursive_distance[0]}-SNP subclusters...")
                result = subprocess.run(f"python3 {script_path} '{args.mat_tree}' '{args.nwk_tree}' -d {next_recursion} -s{samples_in_cluster_str} -v {args_dot_type} -neo -nm -nl -sf {number_part+1} {subsequent_recursions}", shell=True)
                number_start = result.returncode
            logging.debug('%s done recursing, number_start was %s, number_part is %s, number_start now %s', cluster_name, number_start_old, number_part, number_start)
        
        # add in the unclustered samples (outside for loop to avoid writing multiple times)
        # however, don't add to the UUID list, or else persistent cluster IDs will break
        if not args.nolonely:
            lonely = sorted(list(lonely))
            for george in lonely: # W0621, https://en.wikipedia.org/wiki/Lonesome_George
                sample_cluster.append(f"{george}\tlonely\n")
            with open(f"{prefix}_lonely.txt", "a", encoding="utf-8") as unclustered_samples_list: # TODO: again, should this be append or write?
                unclustered_samples_list.writelines(line + '\n' for line in lonely)
            unclustered_as_str = ','.join(lonely)
            cluster_samples.append(f"lonely\t{unclustered_as_str}\n")
            handle_subprocess("Also extracting a tree for lonely samples...",
                f'matUtils extract -i "{args.mat_tree}" -t "LONELY" -s {prefix}_lonely.txt -N {minimum_tree_size}')
            logging.debug(os.listdir('.'))

        if not args.noextraouts:
            logging.info("Writing auspice-style TSV for annotation of clusters...")
            with open(f"{prefix}_cluster_annotation.tsv", "a", encoding="utf-8") as samples_for_annotation:
                samples_for_annotation.writelines(sample_cluster)

            logging.info("Writing auspice-style TSV with cluster UUIDs instead of full names; used for persistent cluster IDs")
            # this one should never include unclustered samples
            with open(f"{prefix}_cluster_UUIDs.tsv", "a", encoding="utf-8") as samples_by_cluster_UUID:
                samples_by_cluster_UUID.writelines(sample_clusterUUID)
            
            logging.info("Writing usher-style TSV for subtree extraction... (although this isn't really the one we will use here)")
            # note that because we are recursing, samples can have more than one subtree assignment, so this can't be fed into usher directly
            with open(f"{prefix}_cluster_extraction.tsv", "a", encoding="utf-8") as clusters_for_subtrees:
                cluster_samples.append("\n") # to avoid skipping last line when read
                clusters_for_subtrees.writelines(cluster_samples)

            logging.info("Writing other little doodads...")
            with open("n_big_clusters", "w", encoding="utf-8") as n_cluster: n_cluster.write(str(n_clusters))
            with open("n_samples_in_clusters", "w", encoding="utf-8") as n_cluded: n_cluded.write(str(n_samples_in_clusters))
            with open("n_samples_processed", "w", encoding="utf-8") as n_processed: n_processed.write(str(total_samples_processed))
            with open("n_unclustered", "w", encoding="utf-8") as n_lonely: n_lonely.write(str(len(lonely)))
    else:
        logging.debug("Not clustering, exiting...")
    
    if not args.nocluster:
        last_cluster_number = number_start
        logging.debug('returning %s @ %sSNP', last_cluster_number, args.distance)
    else:
        last_cluster_number = 99999999999999999999999999999999999999999999999999999999999  # should never get read, but if it does, we'll certainly notice!
    exit(last_cluster_number)


def handle_subprocess(explainer, system_call_as_string):
    logging.info(explainer)
    time.sleep(1)
    logging.debug(system_call_as_string) # pylint: disable=W1203
    time.sleep(1)
    os.system(system_call_as_string)
    time.sleep(1)
    return

def path_to_root(ete_tree, node_name):
    # Browse the tree from a specific leaf to the root
    #logging.debug("Getting path for %s in %s", node_name, type(ete_tree))
    node = ete_tree.search_nodes(name=node_name)[0]
    #logging.debug("Node as found in ete tree: %s", node)
    path = [node]
    while node:
        node = node.up
        path.append(node)
    #logging.debug("path for %s: %s", node_name, path)
    return path

def dist_matrix(tree_to_matrix, samples, args):
    samp_ancs = {}
    #samp_dist = {}
    neighbors = []
    unclustered = set()
    
    #for each input sample, find path to root and branch lengths
    for sample in samples:
    #for sample in progressbar.tqdm(samples, desc="Finding roots and branch lengths"):
        s_ancs = path_to_root(tree_to_matrix, sample)
        samp_ancs[sample] = s_ancs
    
    #create matrix for samples
    matrix = np.full((len(samples),len(samples)), -1)

    for i in range(len(samples)):
    #for i in progressbar.trange(len(samples), desc="Creating matrix"): # trange is a tqdm optimized version of range
        this_samp = samples[i]
        definitely_in_a_cluster = False
        #logging.debug("Checking %s", this_samp)

        for j in range(len(samples)):
            that_samp = samples[j]
            #Future goal: add catch to prevent reiteration of already checked pairs
            if that_samp == this_samp: # self-to-self
                matrix[i][j] = '0'
            elif matrix[i][j] == -1: # ie, we haven't calculated this one yet
                #find lca, add up branch lengths 
                this_path = 0
                that_path = 0
                
                for a in samp_ancs[this_samp]:
                    this_path += a.dist
                    if a in samp_ancs[that_samp]:
                        lca = a
                        this_path -= a.dist
                        #logging.debug(f"  found a in samp_ancs[that_samp], setting this_path")
                        break
                
                for a in samp_ancs[that_samp]:
                    that_path += a.dist
                    if a == lca:
                        #logging.debug(f'  a == lca, setting that_path')
                        that_path -= a.dist
                        break
                
                #logging.debug("  sample %s vs other sample %s: this_path %s, that_path %s", this_samp, that_samp, this_path, that_path)
                total_distance = int(this_path + that_path)
                matrix[i][j] = total_distance
                matrix[j][i] = total_distance
                if not args.nocluster and total_distance <= args.distance:
                    #logging.debug("  %s and %s seem to be in a cluster (%s)", this_samp, that_samp, total_distance)
                    neighbors.append(tuple((this_samp, that_samp)))
                    definitely_in_a_cluster = True
        
        # after iterating through all of j, if this sample is not in a cluster, make note of that
        if not args.nocluster and not definitely_in_a_cluster:
            #logging.debug("  %s is either not in a cluster or clustered early", this_samp)
            #logging.debug(matrix[i])
            second_smallest_distance = np.partition(matrix[i], 1)[1] # second smallest, because smallest is self-self at 0
            if second_smallest_distance <= args.distance:
                #logging.debug("  Oops, %s was already clustered! (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                pass
            else:
                #logging.debug("  %s appears to be truly unclustered (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                unclustered.add(this_samp)
    
    # finished iterating, let's see what our clusters look like
    if not args.nocluster:
        true_clusters = []
        first_iter = True
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
    if args.nocluster:
        true_clusters = None
    #logging.debug("Returning:\n\tsamples:\n%s\n\tmatrix:\n%s\n\ttrue_clusters:\n%s\n\tunclustered:\n%s", samples, matrix, true_clusters, unclustered)
    return samples, matrix, true_clusters, unclustered


if __name__ == "__main__":
    main()
