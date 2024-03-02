version 1.0

import "./matutils_and_friends.wdl" as matWDLlib

workflow Convert_to_Newick {
    input {
        File pb_tree
        String outfile_nwk = "converted.nwk"
    }

    call matWDLlib.convert_to_newick as to_newick {
        input:
            input_mat = pb_tree,
            outfile_nwk = outfile_nwk     
    }

    output {
        File tree_nwk = to_newick.newick_tree
    }

}