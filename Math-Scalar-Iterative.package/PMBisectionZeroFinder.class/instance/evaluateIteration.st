operation
evaluateIteration
	| y |
	result := positiveX + negativeX * 0.5.
	y := functionBlock value: result.
	y <= 0 ifTrue: [ negativeX := result ].
	y >= 0 ifTrue: [ positiveX := result ].
	^ self relativePrecision: (positiveX - negativeX) abs 
	