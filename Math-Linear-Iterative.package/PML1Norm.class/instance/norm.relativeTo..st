operation
norm: aVector relativeTo: reference
	^ (aVector inject: 0 into: [ :sum :each | sum + each abs ]) /
		(reference inject: 0.0 into: [ :sum :each | sum + each abs ])