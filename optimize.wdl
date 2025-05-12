version 1.0

import "./matutils_and_friends.wdl" as matWDLlib

workflow matOptimize {
	input {
		File input_mat
		Int max_hours = 1
		Float min_improvement = 0.00000001
	}

	call matWDLlib.matOptimize as matOptimizer {
		input:
			input_mat = input_mat,
			min_improvement = min_improvement,
			max_hours = max_hours
	}

	output {
		File optimized_tree = matOptimizer.optimized_tree
	}
}