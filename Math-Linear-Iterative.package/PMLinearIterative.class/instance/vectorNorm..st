accessing
vectorNorm: anObject
	(anObject respondsTo: #norm:relativeTo:)
		ifFalse: [ self error: 'Vector norm must respond to norm:relativeTo:' ].
	vectorNorm := anObject