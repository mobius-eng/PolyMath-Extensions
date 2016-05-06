operation
norm: aVector
	^ (aVector inject: 0.0 into: [ :sum :each | each squared + sum ]) sqrt