tests
testDiagonalMatrixSolution
	"Tests simple diagoanl system solution"
	| matrix rhsVector initValue cg result solution preconditioner |
	matrix := DhbMatrix rows: #(#(1 0 0) (0 2 0) (0 0 3)).
	rhsVector := #(1 2 3) asDhbVector.
	solution := #(1 1 1) asDhbVector.
	initValue := #(5 5 5) asDhbVector.
	preconditioner := PMIdentityPreconditioner new.
	cg := PMPreconditionedCG l2NormSolverFor: matrix rightHandSide: rhsVector initialValue: initValue.
	cg preconditioner: preconditioner.
	result := cg evaluate - solution.
	self assert: cg hasConverged.
	self assert: (result * result) sqrt <= DhbFloatingPointMachine new defaultNumericalPrecision.
	