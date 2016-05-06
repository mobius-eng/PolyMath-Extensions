operation
norm: aVector relativeTo: reference
	^ ((aVector inject: 0.0 into: [ :sum :each | each squared + sum ]) /
	(reference inject: 0.0 into: [ :sum :each | each squared + sum ])) sqrt