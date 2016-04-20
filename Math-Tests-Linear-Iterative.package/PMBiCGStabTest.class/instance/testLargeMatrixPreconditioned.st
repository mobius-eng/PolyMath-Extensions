tests
testLargeMatrixPreconditioned
	"This test shows that unpreconditioned CG may converge very slowly. In this particular case,
	the method fails to converge in a small number number of iterations."
	| matrix size rhsVector initValue cg result residual |
	size := 1000.
	matrix := DhbMatrix new: size.
	1 to: size do: [ :i |
		1 to: size do: [ :j |
			matrix at: i at: j put: 0.0.
			i = j ifTrue: [ matrix at: i at: j put: 0.5 + i asFloat sqrt ].
			(i - j) abs = 1 ifTrue: [ matrix at: i at: j put: 1.0 ].
			(i - j) = 100 ifTrue: [ matrix at: i at: j put: 1.0 ] ] ].
	rhsVector := DhbVector new: size withAll: 1.0.
	initValue := DhbVector new: size withAll: 1.0.
	cg := PMPreconditionedBiCGStab ewtRMSSolverFor: matrix rightHandSide: rhsVector initialValue: initValue.
	cg maximumIterations: 20.
	cg preconditioner: (PMDiagonalPreconditioner fromMatrix: matrix).
	result := cg evaluate.
	residual := rhsVector - (matrix * result).
	self assert: cg hasConverged.
	self assert: cg residualError <= 1.
	self assert: (PMErrorWeightVectorNorm new defaultOfSize: rhsVector size; reference: rhsVector; norm: residual) <= 1.
	"self assert: (residual * residual) sqrt < 1e-6."
	