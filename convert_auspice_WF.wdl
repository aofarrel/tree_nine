version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Convert_PB_to_Nextstrain_Auspice {
    input {
        File pb_tree
        Array[File]? metadata_files
    }
    # matUtils extract's help describes the metadata files as follows:
    # "Comma delimited names of tsvs or csvs containing sample identifiers in the first column
    # and an arbitrary number of metadata values in separate columns, including a header line
    # in each file."
    
    String basename = basename(pb_tree, ".pb")
    
    call treenine.convert_to_nextstrain_single as convert_to_nextstrain {
        input:
            input_mat = pb_tree,
            metadata_files = metadata_files,
            outfile_nextstrain = basename + ".json"
    }
    
    output {
        File tree_nextstrain = convert_to_nextstrain.nextstrain_singular_tree
    }
}