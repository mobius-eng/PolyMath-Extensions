accessing
linearOperator: anObject
	"Sets linear operator (matrix or a function) of the equation"
	(anObject respondsTo: #value:)
		ifTrue: [ ^ self setFunction: anObject ].
	(anObject respondsTo: #*)
		ifTrue: [ ^ self setFunction: [ :x | anObject * x ] ].
	self error: 'Unknown linear operator object'.