version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Mask_By_Position {
    input {
        File pb_tree
        File mask_tsv
    }

    call treenine.mask as mask {
        input:
            input_mat = pb_tree,
            mask_tsv = mask_tsv
    }

    output {
        File masked_tree = mask.masked_tree
    }
}