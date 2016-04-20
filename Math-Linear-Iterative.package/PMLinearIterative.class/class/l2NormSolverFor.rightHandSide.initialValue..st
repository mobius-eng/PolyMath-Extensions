instance creation
l2NormSolverFor: linOperator rightHandSide: rhsVector initialValue: initValue
	^ self new linearOperator: linOperator;
		rightHandSideVector: rhsVector;
		vectorNorm: PML2Norm new;
		initialValue: initValue;
		desiredPrecision: DhbFloatingPointMachine new defaultNumericalPrecision; yourself.