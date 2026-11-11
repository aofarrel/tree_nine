version 1.0
import "./matutils_and_friends.wdl" as matUtils

workflow convert_to_taxonium {
    input {
        File pb_tree
        String? outfile
    }

    String out = select_first([outfile, basename(pb_tree, ".pb")])

    call matUtils.convert_to_taxonium as convert_to_tax {
        input:
            input_mat = pb_tree,
            outfile_taxonium = out + ".jsonl.gz"
    }

    output {
        File taxonium_tree = convert_to_tax.taxonium_tree
    }

}