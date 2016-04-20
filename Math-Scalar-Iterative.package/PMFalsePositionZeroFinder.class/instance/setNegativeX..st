as yet unclassified
setNegativeX: aNumber
	| y |
	(y := functionBlock value: aNumber) < 0
		ifFalse: [ self error: 'Function is not negative at x = ', aNumber asString ].
	negativeX := aNumber.
	negativeXValue := y