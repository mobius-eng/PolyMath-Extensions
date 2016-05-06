operation
norm: aVector
	^aVector inject: 0 into: [ :sum :each | sum + each abs ]