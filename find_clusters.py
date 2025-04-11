print("FIND CLUSTERS - VERSION 2.1.0")
SCRIPT_PATH = '/scripts/find_clusters.py'

# Notes:
# * This script is called once for the original clusters, and several times for locally-masked clusters.
# * 000000 is a special "cluster" that represents the entire tree. Its cluster distance is UINT32_MAX.

# Not implemented:
# * Context samples -- this causes matUtils extract to extract more than one subtree at a time. There's probably a way around this,
#   but no one's requested this feature so we won't waste time trying to implement it.

# pylint: disable=too-complex,pointless-string-statement,multiple-statements,wrong-import-position,no-else-return,unnecessary-pass,useless-suppression,global-statement,use-dict-literal

import os
import argparse
import logging
from datetime import date
from itertools import chain
import subprocess
from collections import defaultdict
import bte
import numpy as np
from tqdm import tqdm
import pandas as pd # im sick and tired of polars' restrictions on TSV output

np.set_printoptions(linewidth=np.inf, threshold=15)

UINT8_MAX = np.iinfo(np.uint8).max   # UNSIGNED!
UINT16_MAX = np.iinfo(np.uint16).max # UNSIGNED!
UINT32_MAX = np.iinfo(np.uint32).max # UNSIGNED!
MATRIX_INTEGER_MAX = UINT32_MAX      # can be changed by args
CURRENT_UUID = np.int32(-1)          # SIGNED!!!!!!!!!!!
TODAY = date.today().isoformat()
TYPE_PREFIX = ''
INITIAL_PB_PATH, INITIAL_PB_BTE, INITIAL_SAMPS = None, None, None
BIG_DISTANCE_MATRIX = None     # Distance matrix of 000000
ALL_CLUSTERS = []              # List of all Cluster() objects, including 000000
SAMPLES_IN_ANY_CLUSTER = set() # Set of samples in any cluster, excluding 000000
UNCLUSTERED_SAMPLES = set()    # Set of samples that are not in any cluster excluding 000000
SAMPLE_CLUSTER = ['Sample\tCluster\n']   # Nextstrain-style TSV for annotation
CLUSTER_SAMPLES = ['Cluster\tSamples\n'] # matUtils extract-style TSV for subtrees
LATEST_CLUSTERS = ['latest_cluster_id\tcurrent_date\tcluster_distance\tn_samples\tminimum_tree_size\tsample_ids\n'] # Used by persistent ID script, excludes unclustered
LATEST_SAMPLES = ['sample_id\tcluster_distance\tlatest_cluster_id\n']                                               # Used by persistent ID script, excludes unclustered

class UnionFind:
    def __init__(self):
        self.parent = dict()

    def find(self, item):
        # path compression
        if self.parent.setdefault(item, item) != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, a, b):
        self.parent[self.find(a)] = self.find(b)

class Cluster():
    def __init__(self, UUID: int, samples: list, distance: np.uint32, input_pb, subcluster: bool, track_unclustered: bool, writetree: bool):
        self.str_UUID = self.set_str_UUID(UUID)
        assert len(samples) == len(set(samples))
        self.samples = sorted(samples)
        self.get_subclusters = False if distance == 5 else subcluster
        self.track_unclustered = track_unclustered
        if distance > UINT32_MAX:
            raise ValueError("🔚distance is a value greater than the unsigned-uint32 maximum used when generating matrices; cannot continue")
        self.cluster_distance = np.uint32(distance)
        logging.info("[%s] Hello, I have %s samples: %s", self.debug_name(), len(self.samples), self.samples)
        self.update_globals()
        
        # initalize other stuff
        self.subclusters = []
        self.unclustered = set()
        self.input_pb = input_pb

        # Currently using a 32-bit signed int matrix in hopes of less aggressive RAM usage
        if MATRIX_INTEGER_MAX == UINT8_MAX:
            self.matrix = np.full((len(samples),len(samples)), 0, dtype=np.uint8)  # UNSIGNED!
        elif MATRIX_INTEGER_MAX == UINT16_MAX:
            self.matrix = np.full((len(samples),len(samples)), 0, dtype=np.uint16) # UNSIGNED!
        else:
            self.matrix = np.full((len(samples),len(samples)), 0, dtype=np.uint32) # SIGNED!

        # Updates self.matrix, self.subclusters, and self.unclustered
        if self.cluster_distance == UINT32_MAX:
            self.subclusters = self.dist_matrix_and_get_subclusters(self.input_pb, 20) # None if not get_subclusters
        elif self.cluster_distance == 20:
            self.subclusters = self.dist_matrix_and_get_subclusters(self.input_pb, 10) # None if not get_subclusters
        elif self.cluster_distance == 10:
            self.subclusters = self.dist_matrix_and_get_subclusters(self.input_pb, 5) # None if not get_subclusters
        else:
            # we already forced self.get_subclusters to false if distance is 5; all we're doing here is setting self.matrix
            self.dist_matrix_and_get_subclusters(self.input_pb, 5)
        
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
        if self.cluster_distance != UINT32_MAX:
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
        i_samples = self.samples  # this was sorted() earlier so it should be sorted in matrix
        j_ghost_index = 0
        neighbors = []

        for i, this_samp in enumerate(tqdm(i_samples, desc="Building matrix")):
            definitely_in_a_cluster = False

            # The pool of samples we allow for j shrinks by one with every iteration of i,
            # in order to prevent calculating distances twice. (We can do this only because
            # our matrix is square and we're starting with two equivalent sorted lists.)
            #
            # We keep track of how many shrinks we have done using "j_ghost_index".
            #
            # on first iteration:
            # j_ghost = 1
            # i = [A*, B, C, D]
            # j = [B, C, D]
            # ---> A:B, A:C, and A:D
            # next:
            # j_ghost = 2
            # i = [A, B*, C, D]
            # j = [C, D]
            # ---> B:C and B:D (we already had B:A from previous iteration)
            # next:
            # j_ghost = 3
            # i = [A, B, C*, D]
            # j = [D]
            # ---> C:D (already have C:A and C:B)
            # next:
            # j_ghost = 4
            # i = [A, B, C, D*]
            # j = []
            # ---> finished
            #
            # We need to keep track of how many shrinks we have done with "j_ghost_index"
            # in order to place these calculated values in the correct place on the matrix.
            # j_ghost_index + enumerate(j_samples) = correct index for the j bit of the matrix
            j_ghost_index += 1
            j_samples = self.samples[j_ghost_index:]

            for j, that_samp in enumerate(j_samples):

                j_matrix = j + j_ghost_index
                logging.debug("j %s, j_ghost_index %s, j_matrix %s, that_samp %s", j, j_ghost_index, j_matrix, that_samp)

                this_node, that_node = this_samp, that_samp
                LCA = tree_to_matrix.LCA([this_samp, that_samp])     # type str
                logging.debug("%s:%s LCA is %s", this_samp, that_samp, LCA)
                total_distance = self.sum_paths_to_LCA_plus_overflow_check(tree_to_matrix, this_node, that_node, LCA)
                self.matrix[i][j_matrix], self.matrix[j_matrix][i] = total_distance, total_distance
                if self.get_subclusters and total_distance <= subcluster_distance:
                    logging.debug("  %s and %s seem to be within a %sSNP-cluster (%s)", this_samp, that_samp, subcluster_distance, total_distance)
                    neighbors.append(tuple((this_samp, that_samp)))
                    definitely_in_a_cluster = True

            # Consider samples A, B, C, D, and E. When i = A, j=B, so we calculate their distance, then assign the result to matrix[A][B]
            # and matrix[B][A]. Then j=C, so we get the distance, assign matrix[A][C] and matrix[C][A], etc...
            # Because the j array is shrinking per iteration of i, that can prevent definitely_in_a_cluster from being triggered if it
            # ought to, which is why we need this bit below.
            if self.get_subclusters and not definitely_in_a_cluster:
                second_smallest_distance = np.partition(self.matrix[i], 1)[1] # second smallest, because smallest is self-self at 0
                if second_smallest_distance <= subcluster_distance:
                    #logging.debug("  Oops, %s was already clustered! (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                    pass
                else:
                    #logging.debug("  %s appears to be truly unclustered (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                    if subcluster_distance in (UINT32_MAX, 20): # pylint: disable=else-if-used
                        UNCLUSTERED_SAMPLES.add(this_samp) # only add to global unclustered if it's not in a 20 SNP cluster

        # finished iterating, let's see what our clusters look like
        logging.info("Here is our matrix")
        logging.info(self.matrix)
        subclusters = self.get_true_clusters(neighbors, self.get_subclusters, subcluster_distance) # None if !get_subclusters
        return subclusters

    def sum_paths_to_LCA_plus_overflow_check(self, tree_to_matrix, this_node, that_node, LCA):
        this_path, that_path = 0,0
        while tree_to_matrix.get_node(this_node).id != LCA:
            this_node = tree_to_matrix.get_node(this_node)   # type MATnode
            this_path += this_node.branch_length             # type float
            this_node = this_node.parent.id                  # type str
        while tree_to_matrix.get_node(that_node).id != LCA:
            that_node = tree_to_matrix.get_node(that_node)   # type MATnode
            that_path += that_node.branch_length             # type float
            that_node = that_node.parent.id                  # type str
        total_distance_i64 = this_path + that_path
        if total_distance_i64 > MATRIX_INTEGER_MAX:
            # this is a debug instead of a warning because it happens so often in the uint8 case
            logging.debug("Total distance between %s and %s is %s, greater than integer maximum; will store as %s", this_node, that_node, total_distance_i64, MATRIX_INTEGER_MAX)
            return self.convert_64int_to_whatever(MATRIX_INTEGER_MAX)
        else:
            return self.convert_64int_to_whatever(total_distance_i64)

    def convert_64int_to_whatever(self, python_int64):
        if MATRIX_INTEGER_MAX == UINT8_MAX:
            return np.uint8(python_int64)      # UNSIGNED!
        elif MATRIX_INTEGER_MAX == UINT16_MAX:
            return np.uint16(python_int64)     # UNSIGNED!
        else:
            return np.uint32(python_int64)      # SIGNED!

    def get_true_clusters(self, neighbors, get_subclusters, subcluster_distance):
        # From neighbors we generated while making distance matrix, define (sub)clusters
        # Every cluster here is is a list of sample IDs
        if get_subclusters:
            logging.info("[%s] Looking for subclusters @ %d", self.debug_name(), subcluster_distance)
            #logging.debug("[%s] Got this list of neighbors: %s", self.debug_name(), neighbors)
            true_clusters, truer_clusters = [], []

            uf = UnionFind()
            for a, b in neighbors:
                uf.union(a, b)
            clusters_dict = defaultdict(set)
            for sample in uf.parent:
                root = uf.find(sample)
                clusters_dict[root].add(sample)
            true_clusters = list(clusters_dict.values())
            logging.debug("[%s] Got these clusters: %s", self.debug_name(), true_clusters)
            
            # Since the big refactor there shouldn't be any overlapping clusters. But, to prevent issues in
            # process_clusters.py, we are try to catch that scenario here, and if caught, remove
            # the smaller version.
            all_samples = list(chain.from_iterable(true_clusters))
            if len(all_samples) != len(set(all_samples)): # TODO: make this an assert later
                logging.warning("[%s] Detected overlapping subclusters.", self.debug_name())
                true_clusters = self.deal_with_subcluster_overlap(true_clusters)

            for cluster in true_clusters:
                logging.debug("[%s] For cluster %s in true_clusters %s", self.debug_name(), cluster, true_clusters)
                if subcluster_distance == UINT32_MAX:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), UINT32_MAX, self.input_pb, True, True, True))
                elif subcluster_distance == 20:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 20, self.input_pb, True, False, True))
                elif subcluster_distance == 10:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 10, self.input_pb, True, False, True))
                else:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 5, self.input_pb, False, False, True))
            return truer_clusters
        else:
            return None

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
        logging.debug("[%s] Wrote distance matrix to %s", self.debug_name(), matrix_out)
        if os.path.getsize(matrix_out) < 52428800 and logging.root.level == logging.DEBUG:
            logging.debug("[%s] It looks like this:", self.debug_name())
            with open(matrix_out, "r", encoding='utf-8') as f:
                print(f.read())
        else:
            logging.debug("[%s] And we're not printing it because it's huge")
        if self.cluster_distance == UINT32_MAX:
            global BIG_DISTANCE_MATRIX
            BIG_DISTANCE_MATRIX = self.matrix

    def deal_with_subcluster_overlap(self, tuples_list):
        logging.debug("[%s] got tuples_list %s of type %s", self.debug_name(), tuples_list, type(tuples_list))
        element_to_tuples, conflicts = defaultdict(set), set()
        for i, tup in enumerate(tuples_list):
            for elem in tup:
                element_to_tuples[elem].add(i)
        logging.debug("[%s] Elements mapped to their tuples: %s", self.debug_name(), element_to_tuples)
        for indices in element_to_tuples.values():
            if len(indices) > 1:
                conflicts.update(indices)

        if conflicts:
            # sort conflicts by size -- we only want the bigger one
            conflicting_tuples = sorted(conflicts, key=lambda idx: (len(tuples_list[idx]), idx))
            logging.debug("[%s] Conflicting tuples: %s", self.debug_name(), conflicting_tuples)
            to_remove, seen_elements = set(), set()
            for idx in conflicting_tuples:
                if any(elem in seen_elements for elem in tuples_list[idx]):
                    to_remove.add(idx)
                else:
                    seen_elements.update(tuples_list[idx])
            tuples_list = [tup for i, tup in enumerate(tuples_list) if i not in to_remove]
        logging.debug("[%s] Returning tuples_list: %s", self.debug_name(), tuples_list)
        return tuples_list

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
    global INITIAL_PB_BTE
    INITIAL_PB_BTE = bte.MATree(INITIAL_PB_PATH)
    global INITIAL_SAMPS
    INITIAL_SAMPS = args.samples.split(',') if args.samples else sorted([leaf.id for leaf in INITIAL_PB_BTE.get_leaves()])
    if args.int8:
        global MATRIX_INTEGER_MAX
        MATRIX_INTEGER_MAX = UINT8_MAX

def next_UUID():
    global CURRENT_UUID
    CURRENT_UUID += 1
    return CURRENT_UUID.copy()

def get_all_20_clusters():
    logging.debug("20 clusters are: %s", [cluster.debug_name() for cluster in ALL_CLUSTERS if cluster.cluster_distance == 20])
    return [cluster for cluster in ALL_CLUSTERS if cluster.cluster_distance == np.uint32(20)]

def setup_clustering(distance):
    # We consider the "whole tree" stuff to be its own cluster that always will exist, which we will kick off like this
    # We will not create ANY actual clusters (20, 10, 5) with this function
    new_cluster = Cluster(next_UUID(), INITIAL_SAMPS, distance, INITIAL_PB_BTE, True, True, True)
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
        handle_subprocess("Geting lonely samples' closest relatives...",
            f'matUtils extract -i "{INITIAL_PB_PATH}" --closest-relatives "lonely_closest_relatives.txt"')
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

def find_neighbors(distance_matrix: np.ndarray, sample_names: list, output_tsv: str, plus_unclustered_focus: bool):
    # Note that this REQUIRES the distance matrix and sample_names list to have the same dimensions and sample order.
    # For this reason, to get an output that only focuses on the unclustered samples (whose closest sample may or may
    # not be a clustered sample, ie, we don't want to just rerun this function on an unclustered-only distance matrix),
    # we just remove rows from the pandas dataframe.
    #
    # Currently unused; it seems buggy...
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
    parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')
    parser.add_argument('-d', '--distance', default=20, type=int, help='max distance between samples to identify as clustered')
    parser.add_argument('-rd', '--recursive-distance', type=lambda x: [int(i) for i in x.strip('"').split(',')], help='after identifying --distance cluster, search for subclusters with these distances')
    parser.add_argument('-t', '--type', choices=['BM', 'NB'], type=str.upper, help='BM=backmasked, NB=not-backmasked; will add BM/NB before prefix')
    parser.add_argument('-cn', '--collection-name', default='unnamed', type=str, help='name of this group of samples (do not include a/b prefix)')
    parser.add_argument('-sf', '--startfrom', default=0, type=int, help='the six-digit int part of cluster UUIDs will begin with the next integer after this one')
    parser.add_argument('-i8', '--int8', action='store_true', help='[untested, not recommended] store distance matrix as 8-bit unsigned integers to save as much memory as possible')
    parser.add_argument('-i16', '--int16', action='store_true', help='[untested] store distance matrix as 16-bit unsigned integers to save memory')
    parser.add_argument('-v', '--verbose', action='store_true', help='enable info logging')
    parser.add_argument('-vv', '--veryverbose', action='store_true', help='enable debug logging')

    # mostly used for recursion
    parser.add_argument('-jmatsu', '--justmatrixandthenshutup', action='store_true', help='just generate a matrix and trees for this cluster then exit')
    args = parser.parse_args()
    initial_setup(args)
    if args.justmatrixandthenshutup:
        Cluster(args.collection_name, INITIAL_SAMPS, args.distance, INITIAL_PB_BTE, False, False, False) # will write dmatrix
    else:
        setup_clustering(UINT32_MAX)
        process_unclustered()
        write_output_files()

if __name__ == "__main__":
    main()
    print("🔚Returning")

