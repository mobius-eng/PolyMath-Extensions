as yet unclassified
testNewton2D2
	| solver diff function |
	solver := PMVectorNewtonZeroFinder new.
	function := [ :v |
		| x1 x2 |
		x1 := v at: 1.
		x2 := v at: 2.
		{	(x1+x2) negated exp negated exp - (x2 * (1 + x1 squared)) .
			 x1 * x2 cos + (x2 * x1 sin) - 0.5 } asDhbVector ].
	solver setFunction: function;
		setJacobianOperator: [ :v |
			| x1 x2 x row1 row2 tmp |
			x1 := v at: 1.
			x2 := v at: 2.
			x := x1 + x2.
			tmp := x negated exp.
			row1 := { tmp negated exp * tmp - (2 * x2 * x1) .
						tmp negated exp * tmp - (1 + x1 squared) }.
			row2 := { x2 cos + (x2 * x1 cos) .
						(x1 * x2 sin) negated + x1 sin }.
			DhbMatrix rows: { row1 . row2 } ];
		initialValue: #(0.0 0.0) asDhbVector;
		linearSolver: (PMLUPPlug new).
	solver evaluate.
	self assert: solver hasConverged.
	diff := #(0.3532 0.6061) asDhbVector - solver result.
	self assert: (PML2Norm new norm: diff) <= 1.0e-4.
	self assert: (PMErrorWeightVectorNorm new defaultOfSize: 2; norm: (function value: solver result)) <= 1.