tests
testLargeMatrix
	"This test shows that unpreconditioned CG may converge very slowly. In this particular case,
	the method fails to converge in a small number number of iterations."
	| matrix size rhsVector initValue cg result |
	size := 1000.
	matrix := DhbMatrix new: size.
	1 to: size do: [ :i |
		1 to: size do: [ :j |
			matrix at: i at: j put: 0.
			i = j ifTrue: [ matrix at: i at: j put: 0.5 + i asFloat sqrt ].
			(i - j) abs = 1 ifTrue: [ matrix at: i at: j put: 1 ].
			(i - j) abs = 100 ifTrue: [ matrix at: i at: j put: 1 ] ] ].
	rhsVector := DhbVector new: size withAll: 1.
	initValue := DhbVector new: size withAll: 1.
	cg := PMConjugateGradientL2 new.
	cg linearOperator: matrix; rightHandSideVector: rhsVector; initialValue: initValue.
	cg desiredPrecision: DhbFloatingPointMachine new defaultNumericalPrecision.
	cg maximumIterations: 50.
	result := cg evaluate.
	self assert: cg hasConverged not.
	