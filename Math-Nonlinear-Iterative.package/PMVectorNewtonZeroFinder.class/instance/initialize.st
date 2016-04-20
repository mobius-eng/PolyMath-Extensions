initialization
initialize
	"Initialize lineSearchFunction and lineSearch: it depends on other variables only indirectly,
	can be initialized once at the begining"
	| lineSearchFunctionBlock |
	super initialize.
	desiredPrecision := 1.0.
	lineSearchFunctionBlock := [ :t |
		newFunctionValue := functionBlock value: t * searchStep + result.
		0.5 * (newFunctionValue * newFunctionValue) ].
	lineSearch := PMLineSearchFunctionMinimizer new.
	lineSearch setFunction:  lineSearchFunctionBlock.
	^ self