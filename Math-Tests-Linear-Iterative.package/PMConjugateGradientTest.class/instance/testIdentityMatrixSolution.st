tests
testIdentityMatrixSolution
	"On identity matrix I*x=b, Conjugate Gardient should result in x = b"
	| matrix rhsVector initValue cg result |
	matrix := DhbMatrix rows: #(#(1 0 0) (0 1 0) (0 0 1)).
	rhsVector := #(1 2 3) asDhbVector.
	initValue := #(5 5 5) asDhbVector.
	cg := PMConjugateGradientL2 new.
	cg linearOperator: matrix; rightHandSideVector: rhsVector; initialValue: initValue.
	cg desiredPrecision: DhbFloatingPointMachine new defaultNumericalPrecision.
	result := cg evaluate - rhsVector.
	self assert: cg hasConverged.
	self assert: (result * result) sqrt <= DhbFloatingPointMachine new defaultNumericalPrecision.
	