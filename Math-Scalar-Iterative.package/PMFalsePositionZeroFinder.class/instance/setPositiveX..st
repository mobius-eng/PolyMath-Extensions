as yet unclassified
setPositiveX: aNumber
	| y |
	(y := functionBlock value: aNumber) > 0
		ifFalse: [ self error: 'Function is not positive at x = ', aNumber asString ].
	positiveX := aNumber.
	positiveXValue := y