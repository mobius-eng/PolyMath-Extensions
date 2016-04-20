tests
testSymmetricMatrixSolution
	"Tests a general symmetric matrix"
	| matrix rhsVector initValue cg result solution |
	matrix := DhbMatrix rows: #(#(-7  0  0 -1  3)
										#( 0  5  2  0  0)
										#( 0  2 -3  1  0)
										#(-1  0  1  2 -1)
										#( 3  0  0 -1  1)).
	rhsVector := #(-20 2 0 6 2) asDhbVector.
	solution := #(2 0 1 3 -1) asDhbVector.
	initValue := #(1 1 1 1 1) asDhbVector.
	cg := PMConjugateGradient l2NormSolverFor: matrix rightHandSide: rhsVector initialValue: initValue.
	result := cg evaluate - solution.
	self assert: cg hasConverged.
	self assert: (result * result) sqrt <= DhbFloatingPointMachine new defaultNumericalPrecision.
	