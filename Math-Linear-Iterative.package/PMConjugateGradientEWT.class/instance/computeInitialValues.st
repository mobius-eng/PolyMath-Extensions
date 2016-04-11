operation
computeInitialValues
	super computeInitialValues.
	reference isNil ifTrue: [ reference := rightHandSideVector ].
	absTolerance isNil ifTrue: [ absTolerance := DhbVector new: result size withAll: DhbFloatingPointMachine new defaultNumericalPrecision ].
	relTolerance isNil ifTrue: [ relTolerance := DhbVector new: result size withAll: DhbFloatingPointMachine new defaultNumericalPrecision ].