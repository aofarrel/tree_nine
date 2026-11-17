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
	}

	call matWDLlib.matOptimize as matOptimizer {
		input:
			input_mat = input_mat,
			max_hours = max_hours,
			max_iterations = max_iterations,
			min_improvement = min_improvement,
			node_proportion = node_proportion,
			radius = radius
	}

	String matOptimize_arg_1 = "--max_hours" + max_hours
	String matOptimize_arg_2 = "--max_iterations" + max_iterations
	String matOptimize_arg_3 = "--min_improvement" + min_improvement
	String matOptimize_arg_4 = "--node_proportion" + node_proportion
	String matOptimize_arg_5 = "--radius" + radius

	output {
		File optimized_tree = matOptimizer.optimized_tree
		String matOptimize_args = matOptimize_arg_1+matOptimize_arg_2+matOptimize_arg_3+matOptimize_arg_4+matOptimize_arg_5
	}
}