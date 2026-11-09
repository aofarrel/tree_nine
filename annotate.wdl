version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Annotate {
    input {
        File pb_tree
        File metadata_tsv
        Boolean is_mutation_paths = false
    }
    
    String basename = basename(pb_tree, ".pb")

    if (is_mutation_paths) {
        call treenine.annotate_P as annotate_P {
        input:
            input_mat = pb_tree,
            metadata_tsv = metadata_tsv,
            outfile_mat = basename+"_annotated.pb"
        }
    }

    if (!is_mutation_paths) {
        call treenine.annotate_c as annotate_c {
        input:
            input_mat = pb_tree,
            metadata_tsv = metadata_tsv,
            outfile_mat = basename+"_annotated.pb"
        }
    }
    
    
    
    output {
        File annotated_tree = select_first([annotate_P.annotated_tree, annotate_c.annotated_tree])
    }
}