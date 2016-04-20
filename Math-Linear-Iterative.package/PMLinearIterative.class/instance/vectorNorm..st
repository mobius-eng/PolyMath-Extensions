as yet unclassified
vectorNorm: anObject
	(anObject respondsTo: #norm:) ifFalse: [ self error: 'Vector norm must respond to norm:' ].
	vectorNorm := anObject