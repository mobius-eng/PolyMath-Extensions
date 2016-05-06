tests
testNewton2D
	| solver diff |
	solver := PMVectorNewtonZeroFinder new.
	solver setFunction: [ :v |
		| x1 x2 |
		x1 := v at: 1.
		x2 := v at: 2.
		{	2 * x1 - x2 - x1 negated exp .
			x1 negated + (2 * x2) - x2 negated exp } asDhbVector ];
		setJacobianOperator: [ :v |
			| x1 x2 row1 row2  |
			x1 := v at: 1.
			x2 := v at: 2.
			row1 := { 2 + x1 negated exp . -1 }.
			row2 := { -1 . 2 + x2 negated exp }.
			DhbMatrix rows: { row1 . row2 } ];
		initialValue: #(0.0 0.0) asDhbVector;
		linearSolver: (PMLUPPlug new).
	solver evaluate.
	self assert: solver hasConverged.
	diff := #(0.5671 0.5671) asDhbVector - solver result.
	self assert: (PML2Norm new norm: diff) <= 1.0e-4