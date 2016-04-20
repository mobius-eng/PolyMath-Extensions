as yet unclassified
testDiagonalMatrixSolutionL1
	"Tests simple diagoanl system solution"
	| matrix rhsVector initValue cg result solution |
	matrix := DhbMatrix rows: #(#(1 0 0) (0 2 0) (0 0 3)).
	rhsVector := #(1 2 3) asDhbVector.
	solution := #(1 1 1) asDhbVector.
	initValue := #(5 5 5) asDhbVector.
	cg := PMConjugateGradient l1NormSolverFor: matrix rightHandSide: rhsVector initialValue: initValue.
	result := cg evaluate - solution.
	self assert: cg hasConverged.
	self assert: (PML1Norm new norm: result) <= DhbFloatingPointMachine new defaultNumericalPrecision.
	