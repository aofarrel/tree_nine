version 1.0

import "./matutils_and_friends.wdl" as matWDLlib

workflow matOptimize {
	input {
		File  input_mat
		Int   max_hours        = 1
		Int   max_iterations   = 10000
		Float min_improvement  = 0.00000001
		Int   node_proportion  = 2
		Int   radius           = -1

		Boolean put_args_in_filename = false
	}

	String matOptimize_arg_1 = "__--max_hours__" + max_hours
	String matOptimize_arg_2 = "__--max_iterations__" + max_iterations
	String matOptimize_arg_3 = "__--min_improvement__" + min_improvement
	String matOptimize_arg_4 = "__--node_proportion__" + node_proportion
	String matOptimize_arg_5 = "__--radius__" + radius
	String matOptimize_args = matOptimize_arg_1+matOptimize_arg_2+matOptimize_arg_3+matOptimize_arg_4+matOptimize_arg_5

	String normal_outfile = basename(input_mat, ".pb") + "_optimized.pb"
	String outfile_if_args = basename(input_mat, ".pb") + "__matOptimize" + matOptimize_args + ".pb"

	call matWDLlib.matOptimize as matOptimizer {
		input:
			input_mat = input_mat,
			output_mat = if put_args_in_filename then outfile_if_args else normal_outfile,
			max_hours = max_hours,
			max_iterations = max_iterations,
			min_improvement = min_improvement,
			node_proportion = node_proportion,
			radius = radius
	}

	output {
		File optimized_tree = matOptimizer.optimized_tree
		String utilized_matOptimize_args = matOptimize_args
	}
}