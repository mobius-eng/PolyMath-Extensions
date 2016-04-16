operation
computeInitialValues
	"Private - If no derivative has been defined, take an ad-hoc definition.
	 If no initial value has been defined, take 0 as the starting point (for lack of anything better)."
	result isNil
		ifTrue: [ result := 0].
	derivativeBlock isNil
		ifTrue: [ self error: 'Derivative must be provided. If derivative is unknown, use other methods' ].
	newFunctionValue := functionBlock value: result.