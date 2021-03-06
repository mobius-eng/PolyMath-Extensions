operation
evaluateIteration
	"Compute one step of Newton's zero finding method. Answers the estimated precision."
	| g0 dg0 coefficient |
	"Calculate new Newton step"
	searchStep := (( lastFunctionValue := newFunctionValue) / ( derivativeBlock value: result)) negated.
	"Initialize variables for line search"
	g0 := 0.5 * lastFunctionValue squared.
	dg0 := -2 * lastFunctionValue.
	"Get coefficient of the step from line search: 0 < coefficient <= 1; newFunctionValue is updated as well"
	coefficient := lineSearch setValueAtZero: g0 derivativeAtZero: dg0; evaluate.
	lineSearch hasConverged ifFalse: [ PMStopIterations new signal ].
	searchStep := coefficient * searchStep.
	result := result + searchStep.
	^ self relativePrecision: (searchStep abs max: newFunctionValue abs)
	