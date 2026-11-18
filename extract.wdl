version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Extract {
    input {
        File pb_tree
        Boolean collapse_tree = false

        Boolean subtree = true
        File? samples_in_subtree
        Int? nearest_k
    }
    
    call treenine.extract as extract {
        input:
            input_mat = pb_tree,
            subtree = subtree,
            samples = samples_in_subtree,
            nearest_k = nearest_k,
            collapse_tree = collapse_tree
    }
    
    output {
        File output_tree = extract.output_tree
    }
}