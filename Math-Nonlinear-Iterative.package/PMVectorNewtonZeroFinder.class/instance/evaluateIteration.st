operation
evaluateIteration
	"Compute one step of Newton's zero finding method. Answers the estimated precision."
	| g0 dg0 coefficient |
	"Calculate new Newton step"
	lastFunctionValue := newFunctionValue.
	linearSolver linearOperator: (jacobianOperator value: result);
		rightHandSideVector: lastFunctionValue; initialValue: result.
	searchStep := linearSolver evaluate negated.
	linearSolver hasConverged ifFalse: [ PMStopIterations new signal ].
	"Initialize variables for line search"
	g0 := 0.5 * (lastFunctionValue * lastFunctionValue).
	dg0 := -2 * g0 squared.
	"Get coefficient of the step from line search: 0 < coefficient <= 1; newFunctionValue is updated as well"
	coefficient := lineSearch setValueAtZero: g0 derivativeAtZero: dg0; evaluate.
	lineSearch hasConverged ifFalse: [ PMStopIterations new signal ].
	result := coefficient * searchStep + result.
	^ searchStep rootMeanSquareNormWith: errorWeightVector
	