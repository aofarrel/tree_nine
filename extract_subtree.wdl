version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Extract {
    input {
        File pb_tree
        File samples_in_subtree
        Int? nearest_k
    }
    
    call treenine.extract as extract {
        input:
            input_mat = pb_tree,
            samples = samples_in_subtree,
            nearest_k = nearest_k
    }
    
    output {
        File subtree = extract.subtree
    }
}