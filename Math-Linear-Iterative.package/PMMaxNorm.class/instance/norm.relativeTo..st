operation
norm: aVector relativeTo: reference
	^ (aVector inject: 0 into: [ :sum :each |
		| x |
		x := each abs.
		x > sum ifTrue: [ x ] ifFalse: [ sum ] ]) /
	(reference inject: 0 into: [ :sum :each |
		| x |
		x := each abs.
		x > sum ifTrue: [ x ] ifFalse: [ sum ] ])