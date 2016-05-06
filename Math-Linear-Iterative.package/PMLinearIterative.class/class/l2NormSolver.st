instance creation
l2NormSolver
	^ self new
		vectorNorm: PML2Norm new;
		desiredPrecision: DhbFloatingPointMachine new defaultNumericalPrecision;
		yourself.