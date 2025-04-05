import os
from datetime import date
from collections import defaultdict
import subprocess
import logging
import numpy as np
import bte
from tqdm import tqdm

np.set_printoptions(linewidth=np.inf, threshold=15)
logging.getLogger().setLevel(logging.INFO)
TREE = bte.MATree("/dvol/70.pb")


INT8_MAX = np.iinfo(np.int8).max
INT16_MAX = np.iinfo(np.int16).max
INT32_MAX = np.iinfo(np.int32).max
INT64_MAX = np.iinfo(np.int64).max
MATRIX_INTEGER_MAX = INT32_MAX    # can be changed by args
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


def setup_clustering(distance):
    # We consider the "whole tree" stuff to be its own cluster that always will exist, which we will kick off like this
    # We will not create ANY actual clusters (20, 10, 5) with this function
    INITIAL_SAMPS=TREE.get_leaves()
    INITIAL_SAMPS = [SAMP.id for SAMP in INITIAL_SAMPS]
    new_cluster = Cluster(next_UUID(), INITIAL_SAMPS, distance, True, True, True)
    ALL_CLUSTERS.append(new_cluster)

def handle_subprocess(explainer, system_call_as_string):
    # Wrapper function matUtils subprocesses
    logging.info(explainer)
    logging.debug(system_call_as_string)
    subprocess.run(system_call_as_string, shell=True, check=True)

def next_UUID():
    global CURRENT_UUID
    CURRENT_UUID += 1
    return CURRENT_UUID.copy()

class Cluster():
    def __init__(self, UUID: int, samples: list, distance: np.int32, subcluster: bool, track_unclustered: bool, writetree: bool):
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

        # Currently using a 32-bit signed int matrix in hopes of less aggressive RAM usage
        # I did try switching to unsigned ints by initializing with None but currently that's not worth the trouble
        if MATRIX_INTEGER_MAX == INT8_MAX:
            self.matrix = np.full((len(samples),len(samples)), 0, dtype=np.int8)
        else:
            self.matrix = np.full((len(samples),len(samples)), 0, dtype=np.int32) # TODO: can we make this np.int8 too?

        # Updates self.matrix, self.subclusters, and self.unclustered
        if self.cluster_distance == INT32_MAX:
            self.subclusters = self.new_dist_matrix_and_subcluster(TREE, 20) # None if not get_subclusters
        elif self.cluster_distance == 20:
            self.subclusters = self.new_dist_matrix_and_subcluster(TREE, 10) # None if not get_subclusters
        elif self.cluster_distance == 10:
            self.subclusters = self.new_dist_matrix_and_subcluster(TREE, 5) # None if not get_subclusters
        else:
            # we already forced self.get_subclusters to false if distance is 5; all we're doing here is setting self.matrix
            self.new_dist_matrix_and_subcluster(TREE, 5)
        
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

    def new_dist_matrix_and_subcluster(self, tree_to_matrix, subcluster_distance):
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
                #if self.matrix[i][j] == -1: # should always be the case
                #if that_samp == this_samp: # self-to-self, should never happen
                #    self.matrix[i][j] = np.int8(0)
                #else:
                # find lca, add up branch lengths

                j_matrix = j + j_ghost_index
                logging.debug(f"j {j}, j_ghost_index {j_ghost_index}, j_matrix {j_matrix}, that_samp {that_samp}")

                this_path, that_path = 0, 0
                this_node, that_node = this_samp, that_samp
                LCA = tree_to_matrix.LCA([this_samp, that_samp])     # type str
                logging.debug(f"{this_samp}:{that_samp} LCA is {LCA}")
                while tree_to_matrix.get_node(this_node).id != LCA:
                    this_node = tree_to_matrix.get_node(this_node)   # type MATnode
                    this_path += this_node.branch_length             # type float
                    this_node = this_node.parent.id                  # type str
                while tree_to_matrix.get_node(that_node).id != LCA:
                    that_node = tree_to_matrix.get_node(that_node)   # type MATnode
                    that_path += that_node.branch_length             # type float
                    that_node = that_node.parent.id                  # type str
                logging.debug("  sample %s vs other sample %s: this_path %s, that_path %s", this_samp, that_samp, this_path, that_path)
                if MATRIX_INTEGER_MAX == INT8_MAX:
                    total_distance = self.matrix_i8(this_path, that_path)
                else:
                    total_distance = self.matrix_i32(this_path, that_path, this_samp, that_samp)
                self.matrix[i][j_matrix], self.matrix[j_matrix][i] = total_distance, total_distance

                logging.debug("ASSINGED AS SUCH")
                logging.debug(self.matrix[i][j_matrix])
                logging.debug(self.matrix[j_matrix][i])

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
                    if subcluster_distance in (INT32_MAX, 20): # pylint: disable=else-if-used
                        UNCLUSTERED_SAMPLES.add(this_samp) # only add to global unclustered if it's not in a 20 SNP cluster

        # finished iterating, let's see what our clusters look like
        logging.info(self.matrix)
        subclusters = self.get_true_clusters(neighbors, self.get_subclusters, subcluster_distance) # None if !get_subclusters
        return subclusters


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
                        self.matrix[i][j] = np.int8(0)
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
                        if MATRIX_INTEGER_MAX == INT8_MAX:
                            total_distance = self.matrix_i8(this_path, that_path)
                        else:
                            total_distance = self.matrix_i32(this_path, that_path, this_samp, that_samp)
                        self.matrix[i][j], self.matrix[j][i] = total_distance, total_distance
                        if self.get_subclusters and total_distance <= subcluster_distance:
                            #logging.debug("  %s and %s seem to be within a %sSNP-cluster (%s)", this_samp, that_samp, subcluster_distance, total_distance)
                            neighbors.append(tuple((this_samp, that_samp)))
                            definitely_in_a_cluster = True

            # Consider samples A, B, C, D, and E. When i = A, j=B, so we calculate their distance, then assign the result to matrix[A][B]
            # and matrix[B][A]. Then j=C, so we get the distance, assign matrix[A][C] and matrix[C][A], etc...
            # When i=B and j=A, then we end up at matrix[B][A], which we already filled in with a non -1 value. So we skip it. But that
            # prevents the definitely_in_a_cluster flag from being triggered if it ought to, which is why we need this bit below.
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
        return subclusters

    def matrix_i8(self, this_path, that_path):
        total_distance = this_path + that_path
        if total_distance > MATRIX_INTEGER_MAX:
            # no warning because this will fire constantly and slow things down
            return np.int8(MATRIX_INTEGER_MAX)
        else:
            return np.int8(total_distance)
    
    def matrix_i32(self, this_path, that_path, this_samp, that_samp):
        # Checks for integer overflow when generating distance matrix
        # The expense of printing this warning again and again isn't worth it if you're tryna squeeze performance using
        # unusual integer sizes, so this will only warn on the default case.
        total_distance = this_path + that_path
        if total_distance > MATRIX_INTEGER_MAX:
                logging.warning("Total distance between %s and %s is %s, greater than integer maximum; will store as %s", this_samp, that_samp, total_distance, INT32_MAX)
                return np.int32(MATRIX_INTEGER_MAX)
        else:
            return np.int32(total_distance)

    def get_true_clusters(self, neighbors, get_subclusters, subcluster_distance):
        # From neighbors we generated while making distance matrix, define (sub)clusters
        # Every cluster here is is a list of sample IDs
        if get_subclusters:
            logging.info("[%s] Looking for subclusters @ %d", self.debug_name(), subcluster_distance)
            true_clusters, first_iter, samples_set, samples_list, truer_clusters = [], True, set(), [], []
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
            
            # Ideally, subclusters should have no overlap in samples, but in very large clusters this isn't always true.
            # In theory we can just live with this, but to prevent issues in process_clusters.py, we are try to catch
            # that scenario here. We will remove the smaller of the two subclusters.
            for potential_subcluster in true_clusters:
                samples_set.update(potential_subcluster)
                samples_list.extend(potential_subcluster)
            if len(samples_set) != len(samples_list):
                logging.warning("[%s] Detected overlapping subclusters. Compare set %s vs list %s", self.debug_name(), samples_set, samples_list)
                true_clusters = self.deal_with_subcluster_overlap(true_clusters)

            for cluster in true_clusters:
                logging.info("[%s] For cluster %s in true_clusters %s", self.debug_name(), cluster, true_clusters)
                if subcluster_distance == INT32_MAX:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), INT32_MAX, True, True, True))
                elif subcluster_distance == 20:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 20, True, False, True))
                elif subcluster_distance == 10:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 10, True, False, True))
                else:
                    truer_clusters.append(Cluster(next_UUID(), list(cluster), 5, False, False, True))
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
        logging.debug("[%s] Wrote distance matrix to %s", self.debug_name(), matrix_out)
        if os.path.getsize(matrix_out) < 52428800 and logging.root.level == logging.DEBUG: # 50 MiB
            logging.debug("[%s] It looks like this:")
            with open(matrix_out, "r", encoding='utf-8') as f:
                print(f.read())
        else:
            logging.debug("[%s] And we're not printing it because it's huge")
        if self.cluster_distance == INT32_MAX:
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

setup_clustering(INT32_MAX)
