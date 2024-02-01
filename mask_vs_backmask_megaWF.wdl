version 1.0

import "./mask_vs_backmask_WF.wdl" as MVB_WF
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/main/tasks/processing_tasks.wdl" as ashLib

workflow Big_Cross_Sample_Masking_Test {
    input {
        Array[File] diffs
        File? input_tree
    }
    
    call MVB_WF.Cross_Sample_Masking_Test as one {
        input:
            diffs = diffs,
            input_tree=input_tree,
            out_prefix="01",
            matrix_only_new_samples = true
    }
    
    call MVB_WF.Cross_Sample_Masking_Test as two {
        input:
            diffs = diffs,
            input_tree=input_tree,
            out_prefix="02",
            matrix_only_new_samples = true
    }
    
    call MVB_WF.Cross_Sample_Masking_Test as three {
        input:
            diffs = diffs,
            input_tree=input_tree,
            out_prefix="03",
            matrix_only_new_samples = true
    }
    
    call MVB_WF.Cross_Sample_Masking_Test as four {
        input:
            diffs = diffs,
            input_tree=input_tree,
            out_prefix="04",
            matrix_only_new_samples = true
    }
    
    call ashLib.gather_files as move_files_around {
        input:
            some_files = [
                        one.min_distance_matrix, one.max_taxonium, one.min_taxonium, one.bm_taxonium, one.bm_distance_matrix,
                        two.min_distance_matrix, two.max_taxonium, two.min_taxonium, two.bm_taxonium, two.bm_distance_matrix,
                        three.min_distance_matrix, three.max_taxonium, three.min_taxonium, three.bm_taxonium, three.bm_distance_matrix,
                        four.min_distance_matrix, four.max_taxonium, four.min_taxonium, four.bm_taxonium, four.bm_distance_matrix,
            ]
    }
    
    output {
        Array[File] four_runs = move_files_around.the_same_files
    }
}

