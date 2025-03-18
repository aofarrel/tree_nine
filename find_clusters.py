# FIND CLUSTERS - VERSION 2.0.0
SCRIPT_PATH = '/scripts/find_clusters.py'

# Notes:
# * This script is called once for the original clusters, and several times for locally-masked clusters.
# * 000000 is a special "cluster" that represents the entire tree. Its cluster distance is INT32_MAX.

# Not implemented:
# * Context samples -- this causes matUtils extract to extract more than one subtree at a time. There's probably a way around this,
#   but no one's requested this feature so we won't waste time trying to implement it.

# pylint: disable=too-complex,pointless-string-statement,multiple-statements,wrong-import-position,no-else-return,unnecessary-pass,consider-using-enumerate,useless-suppression,global-statement

import os
import argparse
import logging
from datetime import date
import subprocess
import numpy as np
import ete3
import pandas as pd # im sick and tired of polars' restrictions on TSV output


INT16_MAX = np.iinfo(np.int16).max
INT32_MAX = np.iinfo(np.int32).max
INT64_MAX = np.iinfo(np.int64).max
TODAY = date.today().isoformat()
CURRENT_UUID = np.int32(-1)
TYPE_PREFIX = ''
INITIAL_PB_PATH = None
INITIAL_NWK_PATH = None
INITIAL_NWK_ETE = None
INITIAL_SAMPS = None

BIG_DISTANCE_MATRIX = None     # Distance matrix of 000000
ALL_CLUSTERS = []              # List of all Cluster() objects, including 000000
SAMPLES_IN_ANY_CLUSTER = set() # Set of samples in any cluster, excluding 000000
UNCLUSTERED_SAMPLES = set()    # Set of samples that are not in any cluster excluding 000000
SAMPLE_CLUSTER = ['Sample\tCluster\n']   # Nextstrain-style TSV for annotation
CLUSTER_SAMPLES = ['Cluster\tSamples\n'] # matUtils extract-style TSV for subtrees
LATEST_CLUSTERS = ['latest_cluster_id\tcurrent_date\tcluster_distance\tn_samples\tminimum_tree_size\tsample_ids\n'] # Used by persistent ID script, excludes unclustered
LATEST_SAMPLES = ['sample_id\tcluster_distance\tlatest_cluster_id\n']                                               # Used by persistent ID script, excludes unclustered


class Cluster():
    def __init__(self, UUID: int, samples: list, distance: np.int32, input_nwk, subcluster: bool, track_unclustered: bool, writetree: bool):
        self.str_UUID = self.set_str_UUID(UUID)
        assert len(samples) == len(set(samples))
        self.samples = sorted(samples)
        self.get_subclusters = False if distance == 5 else subcluster
        self.track_unclustered = track_unclustered
        if distance > INT32_MAX:
            raise ValueError("🔚distance is a value greater than the signed-int32 maximum used when generating matrices; cannot continue")
        self.cluster_distance = np.int32(distance)
        logging.info("[%s] Hello, I have %s samples: %s", self.debug_name(), len(self.samples), self.samples)
        self.update_globals()
        
        # initalize other stuff
        self.subclusters = []
        self.unclustered = set()
        # input_nwk is used to calculate this cluster's distance matrix... does it necessarily match the output nwk of this cluster?
        self.input_nwk = input_nwk

        # Currently using a 32-bit signed int matrix in hopes of less aggressive RAM usage
        # I did try switching to unsigned ints by initializing with None but currently that's not worth the trouble
        self.matrix = np.full((len(samples),len(samples)), -1, dtype=np.int32)

        # Updates self.matrix, self.subclusters, and self.unclustered
        if self.cluster_distance == INT32_MAX:
            self.subclusters = self.dist_matrix_and_get_subclusters(self.input_nwk, 20) # None if not get_subclusters
        elif self.cluster_distance == 20:
            self.subclusters = self.dist_matrix_and_get_subclusters(self.input_nwk, 10) # None if not get_subclusters
        elif self.cluster_distance == 10:
            self.subclusters = self.dist_matrix_and_get_subclusters(self.input_nwk, 5) # None if not get_subclusters
        else:
            # we already forced self.get_subclusters to false if distance is 5; all we're doing here is setting self.matrix
            self.dist_matrix_and_get_subclusters(self.input_nwk, 5)
        
        if self.get_subclusters:
            logging.info("[%s] Processed %s samples, found %s subclusters", self.debug_name(), len(self.samples), len(self.subclusters))
        else:
            logging.info("[%s] Processed %s samples (not subclustering)", self.debug_name(), len(self.samples)) # sort of a misnomer since we already clustered...        

        # write distance matrix (and subtree in two formats)
        self.write_dmatrix()
        if writetree:
            self.write_subtrees()

    def set_str_UUID(self, int_UUID):
        return str(int_UUID).zfill(6)

    def update_globals(self):
        # Doesn't set BIG_DISTANCE_MATRIX since we call this function before calling the distance matrix function (and we do that to get
        # some semblance of order, lest the 5SNP clusters end up here first, which would probably be fine I think but a bit weird)
        if self.cluster_distance != INT32_MAX:
            ALL_CLUSTERS.append(self)
            SAMPLES_IN_ANY_CLUSTER.add(sample for sample in self.samples)
            CLUSTER_SAMPLES.append(f"{self.str_UUID}\t{','.join(self.samples)}\n")                         # ⬇️ minimum_tree_size
            LATEST_CLUSTERS.append(f"{self.str_UUID}\t{TODAY}\t{self.cluster_distance}\t{len(self.samples)}\t{len(self.samples)}\t{self.samples}\n")
            for s in self.samples:
                SAMPLE_CLUSTER.append(f"{s}\t{self.str_UUID}\n")
                LATEST_SAMPLES.append(f"{s}\t{self.cluster_distance}\t{self.str_UUID}\n")

    def debug_name(self):
        return f"{self.str_UUID}@{str(self.cluster_distance).zfill(2)}"

    def dist_matrix_and_get_subclusters(self, tree_to_matrix, subcluster_distance):
        # Updates self.matrix, self.subclusters, and self.unclustered
        samples = self.samples # this was sorted() earlier so it should be sorted in matrix
        samp_ancs, neighbors = {}, []

        for sample in samples:
            s_ancs = self.path_to_root(tree_to_matrix, sample)
            samp_ancs[sample] = s_ancs

        for i in range(len(samples)):
            this_samp = samples[i]
            definitely_in_a_cluster = False

            for j in range(len(samples)):
                that_samp = samples[j]
                if self.matrix[i][j] == -1: # ie, we haven't calculated this one yet
                    if that_samp == this_samp: # self-to-self
                        self.matrix[i][j] = np.int32(0)
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

                        #logging.debug("  sample %s vs other sample %s: this_path %s, that_path %s", this_samp, that_samp, this_path, that_path)
                        total_distance = self.matrix_overflow_check(this_path, that_path, this_samp, that_samp)
                        self.matrix[i][j], self.matrix[j][i] = total_distance, total_distance
                        if self.get_subclusters and total_distance <= subcluster_distance:
                            #logging.debug("  %s and %s seem to be within a %sSNP-cluster (%s)", this_samp, that_samp, subcluster_distance, total_distance)
                            neighbors.append(tuple((this_samp, that_samp)))
                            definitely_in_a_cluster = True

            # after iterating through all that_samples (j), if this_sample (i) is not in a cluster, make note of that
            if self.get_subclusters and not definitely_in_a_cluster:
                second_smallest_distance = np.partition(self.matrix[i], 1)[1] # second smallest, because smallest is self-self at 0
                if second_smallest_distance <= subcluster_distance:
                    #logging.debug("  Oops, %s was already clustered! (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                    pass
                else:
                    #logging.debug("  %s appears to be truly unclustered (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                    if subcluster_distance in (INT32_MAX, 20): # pylint: disable=else-if-used
                        UNCLUSTERED_SAMPLES.add(this_samp) # only add to global unclustered if it's not in a 20 SNP cluster

        # finished iterating, let's see what our clusters look like
        subclusters = self.get_true_clusters(neighbors, self.get_subclusters, subcluster_distance) # None if !get_subclusters
        logging.debug("[%s] Matrix:\nmatrix:\n%s", self.debug_name(), self.matrix)
        return subclusters
    
    def matrix_overflow_check(self, this_path, that_path, this_samp, that_samp):
        # Checks for integer overflow when generating distance matrix
        total_distance = this_path + that_path
        if total_distance > INT32_MAX:
            logging.warning("Total distance between %s and %s is %s, greater than signed 32-bit maximum; will store as %s", this_samp, that_samp, total_distance, INT32_MAX)
            return np.int32(INT32_MAX)
        else:
            return np.int32(total_distance)

    def get_true_clusters(self, neighbors, get_subclusters, subcluster_distance):
        # From neighbors we generated while making distance matrix, define (sub)clusters
        # Every cluster here is is a list of sample IDs
        if get_subclusters:
            logging.info("[%s] Looking for subclusters @ %d", self.debug_name(), subcluster_distance)
            true_clusters, first_iter, truer_clusters = [], True, []
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
            for cluster in true_clusters:
                logging.debug("[%s] For cluster %s in true_clusters %s", self.debug_name(), cluster, true_clusters)
                if subcluster_distance == INT32_MAX:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), INT32_MAX, self.input_nwk, True, True, True))
                elif subcluster_distance == 20:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 20, self.input_nwk, True, False, True))
                elif subcluster_distance == 10:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 10, self.input_nwk, True, False, True))
                else:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 5, self.input_nwk, False, False, True))
            return truer_clusters
        else:
            return None

    def path_to_root(self, ete_tree, node_name):
        # Browse the tree from a specific leaf to the root
        node = ete_tree.search_nodes(name=node_name)[0]
        path = [node]
        while node:
            node = node.up
            path.append(node)
        #logging.debug("path for %s: %s", node_name, path)
        return path

    def write_subtrees(self):
        # It would probably more effiecient to extract all subtrees for all clusters at once, rather than one per cluster, but this
        # is easier to implement and keep track of.
        # TODO: also extract JSON version of the tree and add metadata to it (-M metadata_tsv) even though that doesn't go to MR
        tree_outfile = f"{TYPE_PREFIX}{self.str_UUID}" # extension breaks if using -N, see https://github.com/yatisht/usher/issues/389
        assert not os.path.exists(f"{tree_outfile}.nw"), f"Tried to make subtree called {tree_outfile}.nw but it already exists?!"
        with open("temp_extract_these_samps.txt", "w", encoding="utf-8") as temp_extract_these_samps:
            temp_extract_these_samps.writelines(line + '\n' for line in self.samples)
        handle_subprocess(f"Extracting {tree_outfile} pb for {self.str_UUID}...",
            f'matUtils extract -i "{INITIAL_PB_PATH}" -o {tree_outfile}.pb -s temp_extract_these_samps.txt') # DO NOT INCLUDE QUOTES IT BREAKS THINGS
        handle_subprocess(f"Extracting {tree_outfile} nwk for {self.str_UUID}...",
            f'matUtils extract -i "{INITIAL_PB_PATH}" -t {tree_outfile}.nwk -s temp_extract_these_samps.txt') # DO NOT INCLUDE QUOTES IT BREAKS THINGS
        if os.path.exists(f"{tree_outfile}-subtree-1.nw"):
            logging.warning("Generated multiple subtrees for %s, attempting batch rename (this may break things)", self.str_UUID)
            [os.rename(f, f[:-2] + "nwk") for f in os.listdir() if f.endswith(".nw")] # pylint: disable=expression-not-assigned
        else:
            [os.rename(f, f[:-13] + ".nwk") for f in os.listdir() if f.endswith("-subtree-0.nw")] # pylint: disable=expression-not-assigned
        if os.path.exists("subtree-assignments.tsv"):
            os.rename("subtree-assignments.tsv", "lonely-subtree-assignments.tsv")

    def write_dmatrix(self):
        # Write distance matrix. Also update global distance matrix for entire tree if applicable.
        matrix_out = f"{TYPE_PREFIX}{self.str_UUID}_dmtrx.tsv"
        assert not os.path.exists(matrix_out), f"Tried to write {matrix_out} but it already exists?!"
        with open(matrix_out, "a", encoding="utf-8") as outfile:
            outfile.write('sample\t'+'\t'.join(self.samples))
            outfile.write("\n")                # enumerate causes some type issues, just stick with range(len()) for now
            for k in range(len(self.samples)): # pylint: disable=consider-using-enumerate
                line = [str(int(count)) for count in self.matrix[k]]
                outfile.write(f'{self.samples[k]}\t' + '\t'.join(line) + '\n')
        logging.debug("[%s] Wrote distance matrix: %s", self.debug_name(), matrix_out)
        if self.cluster_distance == INT32_MAX:
            global BIG_DISTANCE_MATRIX
            BIG_DISTANCE_MATRIX = self.matrix



########## Global Functions ###########

def initial_setup(args):
    logging.basicConfig(level=logging.DEBUG if args.veryverbose else logging.INFO if args.verbose else logging.WARNING)
    global TYPE_PREFIX
    if args.type == 'BM':
        TYPE_PREFIX = 'b' # for "backmasked"
    elif args.type == 'NB':
        TYPE_PREFIX = 'a' # for... uh... Absolutelynotbackmasked
    else:
        TYPE_PREFIX = ''
    global INITIAL_PB_PATH
    INITIAL_PB_PATH = args.mat_tree
    global INITIAL_NWK_PATH
    INITIAL_NWK_PATH = args.nwk_tree
    global INITIAL_NWK_ETE
    INITIAL_NWK_ETE = ete3.Tree(INITIAL_NWK_PATH, format=1)
    global INITIAL_SAMPS
    INITIAL_SAMPS = args.samples.split(',') if args.samples else sorted([leaf.name for leaf in INITIAL_NWK_ETE])

def next_UUID():
    global CURRENT_UUID
    CURRENT_UUID += 1
    return CURRENT_UUID.copy()

def get_all_20_clusters():
    logging.debug(cluster.debug_name for cluster in ALL_CLUSTERS)
    logging.debug("20 clusters are: %s", [cluster.debug_name() for cluster in ALL_CLUSTERS if cluster.cluster_distance == 20])
    return [cluster for cluster in ALL_CLUSTERS if cluster.cluster_distance == np.int32(20)]

def setup_clustering(distance):
    # We consider the "whole tree" stuff to be its own cluster that always will exist, which we will kick off like this
    # We will not create ANY actual clusters (20, 10, 5) with this function
    new_cluster = Cluster(next_UUID(), INITIAL_SAMPS, distance, INITIAL_NWK_ETE, True, True, True)
    ALL_CLUSTERS.append(new_cluster)

def process_unclustered():
    # Should not be called if justmatrixandthenshutup
    lonely = sorted(list(UNCLUSTERED_SAMPLES))
    for george in sorted(list(lonely)): # W0621, https://en.wikipedia.org/wiki/Lonesome_George
        SAMPLE_CLUSTER.append(f"{george}\tlonely\n")
    with open(f"{TYPE_PREFIX}_lonely.txt", "w", encoding="utf-8") as unclustered_samples_list:
        unclustered_samples_list.writelines(line + '\n' for line in lonely)
    CLUSTER_SAMPLES.append(f"lonely\t{','.join(lonely)}\n")
    if len(lonely) > 0:
        handle_subprocess("Extracting a tree for lonely samples...",
            f'matUtils extract -i "{INITIAL_PB_PATH}" -t "LONELY" -s {TYPE_PREFIX}_lonely.txt -N {len(lonely)}')
        os.rename("subtree-assignments.tsv", "lonely-subtree-assignments.tsv")
        [os.rename(f, f[:-2] + "nwk") for f in os.listdir() if f.endswith(".nw")] # pylint: disable=expression-not-assigned
    else:
        logging.info("Could not find any unclustered samples")

def handle_subprocess(explainer, system_call_as_string):
    # Wrapper function matUtils subprocesses
    logging.info(explainer)
    logging.debug(system_call_as_string)
    subprocess.run(system_call_as_string, shell=True, check=True)

def write_output_files():
    # Previously we used to use CLUSTER_SAMPLES for usher extraction, but since samples can have more than one subtree
    # assignment, we don't do that anymore. We also previously had two SAMPLE_CLUSTER files, one of which was only UUIDs
    # (from back when UUIDs != internal cluster names) and excluded unclustered samples, but we don't have that file
    # anymore either because latest_samples.tsv (which also excludes unclustered samples) is used instead.
    with open("cluster_annotation_workdirIDs.tsv", "a", encoding="utf-8") as samples_for_annotation:
        samples_for_annotation.writelines(SAMPLE_CLUSTER)
    with open("latest_clusters.tsv", "w", encoding="utf-8") as current_clusters:  # TODO: eventually add old/new samp information
        current_clusters.writelines(LATEST_CLUSTERS)
    with open("latest_samples.tsv", "w", encoding="utf-8") as latest_samples:  # TODO: EVENTUALLY ADD OLD/NEW SAMP INFORMATION
        latest_samples.writelines(LATEST_SAMPLES)
    with open("n_big_clusters", "w", encoding="utf-8") as n_cluster: n_cluster.write(str(len(get_all_20_clusters())))
    with open("n_samples_in_clusters", "w", encoding="utf-8") as n_cluded: n_cluded.write(str(len(SAMPLES_IN_ANY_CLUSTER)))
    with open("n_samples_processed", "w", encoding="utf-8") as n_processed: n_processed.write(str(len(INITIAL_SAMPS)))
    with open("n_unclustered", "w", encoding="utf-8") as n_lonely: n_lonely.write(str(len(UNCLUSTERED_SAMPLES)))
    find_neighbors(BIG_DISTANCE_MATRIX, INITIAL_SAMPS, "all_neighbors.tsv", plus_unclustered_focus=True)

def find_neighbors(distance_matrix: np.ndarray, sample_names: list, output_tsv: str, plus_unclustered_focus: bool):
    # Note that this REQUIRES the distance matrix and sample_names list to have the same dimensions and sample order.
    # For this reason, to get an output that only focuses on the unclustered samples (whose closest sample may or may
    # not be a clustered sample, ie, we don't want to just rerun this function on an unclustered-only distance matrix),
    # we just remove rows from the pandas dataframe.
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
    df.to_csv(output_tsv, sep="\t", index=False)
    if plus_unclustered_focus:
        filtered_rows = df[df["sample"].isin(UNCLUSTERED_SAMPLES)]
        filtered_rows.to_csv("unclustered_neighbors.tsv", sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description="Clusterf...inder")
    parser.add_argument('mat_tree', type=str, help='input MAT (.pb)')
    parser.add_argument('nwk_tree', type=str, help='input nwk (.nwk)')
    parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')
    parser.add_argument('-d', '--distance', default=20, type=int, help='max distance between samples to identify as clustered')
    parser.add_argument('-rd', '--recursive-distance', type=lambda x: [int(i) for i in x.strip('"').split(',')], help='after identifying --distance cluster, search for subclusters with these distances')
    parser.add_argument('-t', '--type', choices=['BM', 'NB'], type=str.upper, help='BM=backmasked, NB=not-backmasked; will add BM/NB before prefix')
    parser.add_argument('-cn', '--collection-name', default='unnamed', type=str, help='name of this group of samples (do not include a/b prefix)')
    parser.add_argument('-sf', '--startfrom', default=0, type=int, help='the six-digit int part of cluster UUIDs will begin with the next integer after this one')
    parser.add_argument('-v', '--verbose', action='store_true', help='enable info logging')
    parser.add_argument('-vv', '--veryverbose', action='store_true', help='enable debug logging')

    # mostly used for recursion
    parser.add_argument('-jmatsu', '--justmatrixandthenshutup', action='store_true', help='just generate a matrix and trees for this cluster then exit')
    args = parser.parse_args()
    initial_setup(args)
    if args.justmatrixandthenshutup:
        Cluster(args.collection_name, INITIAL_SAMPS, args.distance, INITIAL_NWK_ETE, False, False, False) # will write dmatrix
    else:
        setup_clustering(INT32_MAX)
        process_unclustered()
        write_output_files()

if __name__ == "__main__":
    main()
    print("🔚Returning")
