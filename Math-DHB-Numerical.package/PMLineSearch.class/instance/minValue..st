accessing
minValue: anObject
	(anObject isNumber and: [ anObject > 0 ])
		ifFalse: [ self error: 'minimum accepted function value must be a positive number' ].
	minValue := anObject