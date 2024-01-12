version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Convert {
    input {
        File pb_tree
        
        Boolean newick = false
        Boolean nextstrain = false
        Boolean taxonium = true
    }
    
    String basename = basename(pb_tree, ".pb")
    
    if(nextstrain) {
        call treenine.convert_to_nextstrain_single as nex {
            input:
                input_mat = pb_tree,
                outfile_nextstrain = basename + ".json"
        }
    }
    
    if(newick) {
        call treenine.convert_to_newick as new {
            input:
                input_mat = pb_tree,
                outfile_nwk = basename + ".nwk"
        }   
    }
    
    if(taxonium) {
        call treenine.convert_to_taxonium as tax {
            input:
                input_mat = pb_tree,
                outfile_taxonium = basename + ".jsonl.gz"
        }
    }
    

    output {
        File? newick_tree = new.newick_tree
        File? nextstrain_tree = nex.nextstrain_singular_tree
        File? taxonium_tree = tax.taxonium_tree
    }
}