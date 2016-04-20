tests
testIdentityMatrixSolution
	"On identity matrix I*x=b, Conjugate Gardient should result in x = b"
	| matrix rhsVector initValue cg result |
	matrix := DhbMatrix rows: #(#(1 0 0) (0 1 0) (0 0 1)).
	rhsVector := #(1 2 3) asDhbVector.
	initValue := #(5 5 5) asDhbVector.
	cg := PMBiCGStab l1NormSolverFor: matrix rightHandSide: rhsVector initialValue: initValue.
	result := cg evaluate - rhsVector.
	self assert: cg hasConverged.
	self assert: (PML1Norm new norm: result) <= DhbFloatingPointMachine new defaultNumericalPrecision.
	