operation
computeInitialValues
	"Private - If no jacobian has been defined, signal an error.
	 If no initial value has been defined, signal an error."
	result isNil
		ifTrue: [ self error: 'Initial value must be provided' ].
	searchStep := DhbVector new: result size withAll: 0.0.
	jacobianOperator isNil
		ifTrue: [ self error: 'Jacobian must be provided. If derivative is unknown, use other methods' ].
	errorWeightVector isNil
		ifTrue: [ errorWeightVector := PMErrorWeightVectorNorm new defaultOfSize: result size ].
	newFunctionValue := functionBlock value: result.