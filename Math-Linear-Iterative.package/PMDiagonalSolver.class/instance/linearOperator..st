as yet unclassified
linearOperator: anObject
	(anObject  respondsTo: #at:)
		ifFalse: [ self error: 'Diagonal matrix must respond to #at: message' ].
	functionBlock := anObject