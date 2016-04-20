as yet unclassified
norm: aVector
	^ aVector inject: 0 into: [ :sum :each |
		| x |
		x := each abs.
		x > sum ifTrue: [ x ] ifFalse: [ sum ] ]