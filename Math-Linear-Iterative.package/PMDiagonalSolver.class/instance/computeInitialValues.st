as yet unclassified
computeInitialValues
	functionBlock isNil ifTrue: [ self error: 'Matrix was not set' ].
	rightHandSideVector isNil ifTrue: [ self error: 'Right hand side vector was not set' ]