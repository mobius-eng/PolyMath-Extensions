operation
evaluate
	"Perform the iteration until either the desired precision is attained or the number of iterations exceeds the maximum."
	iterations := 0.
	^ self privateEvaluate
	