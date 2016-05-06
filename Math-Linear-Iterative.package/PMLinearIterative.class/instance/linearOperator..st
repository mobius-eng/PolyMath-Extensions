accessing
linearOperator: anObject
	"Sets linear operator (matrix or a function) of the equation"
	(anObject respondsTo: #value:)
		ifTrue: [ ^functionBlock := anObject ].
	(anObject respondsTo: #*)
		ifTrue: [ ^functionBlock := [ :x | anObject * x ] ].
	self error: 'Unknown linear operator object'.