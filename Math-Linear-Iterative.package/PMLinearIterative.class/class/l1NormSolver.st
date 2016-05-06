instance creation
l1NormSolver
	^ self new
		vectorNorm: PML1Norm new;
		desiredPrecision: DhbFloatingPointMachine new defaultNumericalPrecision;
		yourself.