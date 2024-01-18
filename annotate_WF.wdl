version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Annotate {
    input {
        File pb_tree
        File metadata_tsv
    }
    
    String basename = basename(pb_tree, ".pb")
    
    call treenine.annotate as annotate {
        input:
            input_mat = pb_tree,
            metadata_tsv = metadata_tsv,
            outfile_annotated = basename+"_annotated.pb"
    }
    
    output {
        File annotated_tree = annotate.annotated_tree
    }
}