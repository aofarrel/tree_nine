version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Mask_By_Subtree_Position {
    input {
        File pb_tree
        File mask_tsv
        File highly_clonal_samples
    }
    
    call treenine.extract as extract {
        input:
            input_mat = pb_tree,
            samples = highly_clonal_samples
    }

    call treenine.mask as mask {
        input:
            input_mat = extract.subtree,
            mask_tsv = mask_tsv
    }

    output {
        File masked_tree = mask.masked_tree
    }
}