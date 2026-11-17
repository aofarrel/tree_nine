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

	output {
		File optimized_tree = matOptimizer.optimized_tree
	}
}