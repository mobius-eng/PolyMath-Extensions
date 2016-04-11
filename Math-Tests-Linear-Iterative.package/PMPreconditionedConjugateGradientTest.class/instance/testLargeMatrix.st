tests
testLargeMatrix
	"Tests simple diagoanl system solution"
	| matrix size rhsVector initValue cg result preconditioner tmp |
	size := 1000.
	matrix := DhbMatrix new: size.
	1 to: size do: [ :i |
		1 to: size do: [ :j |
			matrix at: i at: j put: 0.0.
			i = j ifTrue: [ matrix at: i at: j put: 0.5 + i asFloat sqrt ].
			(i - j) abs = 1 ifTrue: [ matrix at: i at: j put: 1.0 ].
			(i - j) abs = 100 ifTrue: [ matrix at: i at: j put: 1.0 ] ] ].
	rhsVector := DhbVector new: size withAll: 1.0.
	initValue := DhbVector new: size withAll: 1.0.
	preconditioner := PMDiagonalPreconditioner fromMatrix: matrix.
	cg := PMPreconditionedCGL2 new.
	cg linearOperator: matrix; rightHandSideVector: rhsVector; initialValue: initValue.
	cg desiredPrecision: 1e-6.
	cg preconditioner: preconditioner; maximumIterations: 20.
	result := cg evaluate.
	self assert: cg hasConverged.
	tmp := matrix * result - rhsVector.
	self assert: (tmp * tmp) sqrt <= 1e-6.
	