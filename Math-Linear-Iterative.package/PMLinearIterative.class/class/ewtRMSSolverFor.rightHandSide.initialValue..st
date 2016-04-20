instance creation
ewtRMSSolverFor: linOperator rightHandSide: rhsVector initialValue: initValue
	^ self new linearOperator: linOperator;
		rightHandSideVector: rhsVector;
		vectorNorm: ((PMErrorWeightVectorNorm new)
							defaultOfSize: rhsVector size;
							reference: rhsVector; yourself);
		initialValue: initValue;
		desiredPrecision: 1.0; yourself.