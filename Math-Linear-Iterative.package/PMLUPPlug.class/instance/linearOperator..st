accessing
linearOperator: matrix
	lupDecomposition := nil.
	precision := 2.0.
	(matrix respondsTo: #lupDecomposition ) ifTrue: [ lupDecomposition := matrix lupDecomposition ].
	((matrix isKindOf:#ArrayedCollection) and: [ matrix first isKindOf: #ArrayedCollection ])
		ifTrue: [ lupDecomposition := DhbLUPDecomposition equations: matrix ].
	lupDecomposition isNil ifTrue: [ self error: 'Invalid matrix' ]