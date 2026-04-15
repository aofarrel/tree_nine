version 1.0

workflow Test_Distance_Matrix_Methods {
    input {
        File tree_nwk
        File tree_pb
        File samples
        Int memory_gb = 32
        Int cpu = 12
    }

    call create_distance_matrices {
        input:
            tree_nwk = tree_nwk,
            tree_pb = tree_pb,
            samples = samples,
            memory = memory_gb,
            cpu = cpu
    }

    meta {
        description: "Let's see if PhyloDM really is as fast as they say"
    }
}

task create_distance_matrices {
    input {
        File tree_nwk
        File tree_pb
        File samples
        Int memory
        Int cpu
    }

    command <<<
        set -eux pipefail
        pip install phylodm
        
        python3 << CODE
        TREE_AS_NWK = "~{tree_nwk}"
        TREE_AS_PB = "~{tree_pb}"
        MATUTILS_SUMMARY_SAMPLES = "~{samples}"
        UNCLUSTERED_SAMPLES = set()

        import sys
        import csv
        import time
        import logging
        import numpy as np
        import bte
        from phylodm import PhyloDM
        import dendropy

        np.set_printoptions(threshold=sys.maxsize)
        logging.basicConfig(
            format='[%(asctime)s] %(levelname)s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')
        UINT8_MAX = np.iinfo(np.uint8).max   # UNSIGNED!
        UINT16_MAX = np.iinfo(np.uint16).max # UNSIGNED!
        UINT32_MAX = np.iinfo(np.uint32).max # UNSIGNED!
        UINT64_MAX = np.iinfo(np.uint64).max # UNSIGNED!
        MATRIX_INTEGER_MAX = UINT64_MAX


        with open(MATUTILS_SUMMARY_SAMPLES) as f: next(f); SAMPLES = sorted(line.split("\t")[0] for line in f)

        def main():
            # PhyloDM's native load_from_newick_path() seems to hate nwks from UShER, but
            # for some reason it has no issue with them if passed through dendropy first?
            dendro_tree = dendropy.Tree.get_from_path(TREE_AS_NWK, schema='newick')
            phylodm_tree = PhyloDM.load_from_dendropy(dendro_tree)
            matrix_start_time, matrix_name = time.time(), "PHYLODM"
            phylodm_matrix = phylodm_tree.dm(norm=False)
            logging.warning("[%s] Finished calculating matrix samples in %.2f sec", matrix_name, time.time() - matrix_start_time)

            bte_tree = bte.MATree(TREE_AS_PB)
            bte_matrix_simple = dist_matrix_and_get_subclusters(bte_tree, UINT32_MAX, SAMPLES, get_subclusters=False, matrix_name="SIMPLE BTE")
            bte_matrix_checked = dist_matrix_and_get_subclusters(bte_tree, UINT32_MAX, SAMPLES, get_subclusters=True, matrix_name="DOUBLECHECK BTE")

            assert bte_matrix_simple.shape == bte_matrix_checked.shape
            assert bte_matrix_simple.shape == phylodm_matrix.shape
            print(f"Matrixes all have the same shape: {phylodm_matrix.shape}")
            
            compare_matrices(bte_matrix_simple, bte_matrix_checked, "bte_matrix_simple", "bte_matrix_checked", SAMPLES)
            compare_matrices(bte_matrix_checked, bte_matrix_simple, "bte_matrix_checked", "bte_matrix_simple", SAMPLES)

            # Check labels - note that phyloDM inserts spaces where the tree (and BTE) has underscores!
            phylodm_labels = phylodm_tree.taxa()
            bte_labels = SAMPLES
            assert len(phylodm_labels) == len(bte_labels)
            print(f"Matrixes all have the same number of labels: {len(phylodm_labels)}")
            fixed_labels = list()
            for label in phylodm_labels:
                fixed_labels.append(label.replace(" ", "_"))
            phylodm_labels = fixed_labels

            if sorted(phylodm_labels) != bte_labels:
                print(f"Total number of labels: {len(phylodm_labels)}")
                print(f"Number of intersecting labels: {len(set(phylodm_labels).intersection(set(bte_labels)))}")
                print(f"Symmetric difference: {set(phylodm_labels).symmetric_difference(set(bte_labels))}")
                print("PhyloDM - bte")
                print(set(phylodm_labels) - set(bte_labels))
                print("bte - PhyloDM")
                print(set(bte_labels) - set(phylodm_labels))

            # 1. Create a mapping from labels to indices for the second matrix
            label_to_idx_b = {label: i for i, label in enumerate(bte_labels)}

            # 2. Get the new order of indices based on phylodm_labels
            # (Assumes phylodm_labels and bte_labels contain the same set of strings)
            new_indices = [label_to_idx_b[label] for label in phylodm_labels]

            # 3. Reorder Matrix B to match Matrix A
            # np.ix_ creates a mesh so we reorder rows AND columns
            bte_matrix_simple_aligned = bte_matrix_simple[np.ix_(new_indices, new_indices)]

            compare_matrices(phylodm_matrix, bte_matrix_simple_aligned, "phylodm_matrix", "bte_matrix_simple_aligned", SAMPLES)


        def compare_matrices(matrix_a, matrix_b, a_name, b_name, labels):
            diff_norm = np.linalg.norm(matrix_a - matrix_b)
            print(f"[{a_name}, {b_name}] Frobenius Norm: {diff_norm}")
            if diff_norm == 0.0:
                print(f"[{a_name}, {b_name}] appear to be functionally indentical")
                return
            print(labels)
            print(matrix_a)
            print(labels)
            print(matrix_b)
            sample_differences = np.abs(matrix_a - matrix_b).mean(axis=1)
            silliest_samples = sorted(zip(labels, sample_differences), key=lambda x: x[1], reverse=True)
            print(f"[{a_name}, {b_name}] Silliest (most different) samples: {silliest_samples}")


        def dist_matrix_and_get_subclusters(tree_to_matrix, subcluster_distance, samples: list, get_subclusters, matrix_name="SOME_MATRIX"):
            matrix_start_time = time.time()
            i_samples = samples  # this was sorted() earlier so it should be sorted in matrix
            j_ghost_index = 0
            neighbors = []

            # Currently using a 32-bit unsigned int matrix in hopes of less aggressive RAM usage
            if MATRIX_INTEGER_MAX == UINT8_MAX:
                matrix = np.full((len(samples),len(samples)), 0, dtype=np.uint8)  # UNSIGNED!
            elif MATRIX_INTEGER_MAX == UINT16_MAX:
                matrix = np.full((len(samples),len(samples)), 0, dtype=np.uint16) # UNSIGNED!
            else:
                matrix = np.full((len(samples),len(samples)), 0, dtype=np.uint32) # UNSIGNED!

            for i, this_samp in enumerate(i_samples):
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
                j_samples = samples[j_ghost_index:]

                for j, that_samp in enumerate(j_samples):

                    j_matrix = j + j_ghost_index
                    logging.debug("j %s, j_ghost_index %s, j_matrix %s, that_samp %s", j, j_ghost_index, j_matrix, that_samp)

                    this_node, that_node = this_samp, that_samp
                    LCA = tree_to_matrix.LCA([this_samp, that_samp])     # type str
                    logging.debug("%s:%s LCA is %s", this_samp, that_samp, LCA)
                    total_distance = sum_paths_to_LCA_plus_overflow_check(tree_to_matrix, this_node, that_node, LCA)
                    matrix[i][j_matrix], matrix[j_matrix][i] = total_distance, total_distance
                    if get_subclusters and total_distance <= subcluster_distance:
                        logging.debug("  %s and %s seem to be within a %sSNP-cluster (%s)", this_samp, that_samp, subcluster_distance, total_distance)
                        neighbors.append(tuple((this_samp, that_samp)))
                        definitely_in_a_cluster = True

                # Consider samples A, B, C, D, and E. When i = A, j=B, so we calculate their distance, then assign the result to matrix[A][B]
                # and matrix[B][A]. Then j=C, so we get the distance, assign matrix[A][C] and matrix[C][A], etc...
                # Because the j array is shrinking per iteration of i, that can prevent definitely_in_a_cluster from being triggered if it
                # ought to, which is why we need this bit below.
                if get_subclusters and not definitely_in_a_cluster:
                    second_smallest_distance = np.partition(matrix[i], 1)[1] # second smallest, because smallest is self-self at 0
                    if second_smallest_distance <= subcluster_distance:
                        #logging.debug("  Oops, %s was already clustered! (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                        pass
                    else:
                        #logging.debug("  %s appears to be truly unclustered (closest sample is %s SNPs away)", this_samp, second_smallest_distance)
                        if subcluster_distance in (UINT32_MAX, 20): # pylint: disable=else-if-used  # only add to global unclustered if it's not in a 20 SNP cluster
                            if this_samp in INITIAL_SAMPS:
                                UNCLUSTERED_SAMPLES.add(this_samp) # attempt to fix https://github.com/aofarrel/tree_nine/issues/41

            # finished iterating, let's see what our clusters look like
            #logging.info("Here is our matrix")
            #logging.info(self.matrix)
            # This doesn't print len(self.samples) because that was printed earlier already
            logging.warning("[%s] Finished calculating matrix samples in %.2f sec", matrix_name, time.time() - matrix_start_time)
            return matrix
            #subclusters = self.get_true_clusters(neighbors, self.get_subclusters, subcluster_distance) # None if !get_subclusters
            #return subclusters

        def sum_paths_to_LCA_plus_overflow_check(tree_to_matrix, this_node, that_node, LCA):
            this_path, that_path = 0,0
            while tree_to_matrix.get_node(this_node).id != LCA:
                this_node = tree_to_matrix.get_node(this_node)   # type MATnode
                this_path += this_node.branch_length             # type float
                logging.debug(f"{this_node.id}x{this_node.parent.id} branch length: {this_path}")
                this_node = this_node.parent.id                  # type str
            while tree_to_matrix.get_node(that_node).id != LCA:
                that_node = tree_to_matrix.get_node(that_node)   # type MATnode
                that_path += that_node.branch_length             # type float
                logging.debug(f"{that_node.id}x{that_node.parent.id} branch length: {that_path}")
                that_node = that_node.parent.id                  # type str
            total_distance_i64 = this_path + that_path
            logging.debug(f"Total distance is {total_distance_i64}")
            if total_distance_i64 > MATRIX_INTEGER_MAX:
                # this is a debug instead of a warning because it happens so often in the uint8 case
                logging.debug("Total distance between %s and %s is %s, greater than integer maximum; will store as %s", this_node, that_node, total_distance_i64, MATRIX_INTEGER_MAX)
                return convert_64int_to_whatever(MATRIX_INTEGER_MAX)
            else:
                return convert_64int_to_whatever(total_distance_i64)

        def convert_64int_to_whatever(python_int64):
            if MATRIX_INTEGER_MAX == UINT8_MAX:
                return np.uint8(python_int64)      # UNSIGNED!
            elif MATRIX_INTEGER_MAX == UINT16_MAX:
                return np.uint16(python_int64)     # UNSIGNED!
            else:
                return np.uint32(python_int64)     # UNSIGNED!

        if __name__=="__main__":
            main()

        CODE
    >>>

    runtime {
        bootDiskSizeGb: 15
        cpu: cpu
        disks: "local-disk " + 150 + " SSD"
        docker: "ashedpotatoes/usher-plus:0.6.4ash_2"
        memory: memory + " GB"
    }
}