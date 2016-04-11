operation
updateResult: nextX
	" Move x1 -> x2; g(x1) -> g(x2) and assign x1 = nextX; g(x1) = g(nextX) "
	previousResult := result.
	valueAtPreviousResult := valueAtResult.
	result := nextX.
	valueAtResult := functionBlock value: result.	