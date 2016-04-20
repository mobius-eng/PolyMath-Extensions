instance creation
l1NormSolverFor: linOperator rightHandSide: rhsVector initialValue: initValue
	^ self new linearOperator: linOperator;
		rightHandSideVector: rhsVector;
		vectorNorm: PML1Norm new;
		initialValue: initValue;
		desiredPrecision: DhbFloatingPointMachine new defaultNumericalPrecision; yourself.