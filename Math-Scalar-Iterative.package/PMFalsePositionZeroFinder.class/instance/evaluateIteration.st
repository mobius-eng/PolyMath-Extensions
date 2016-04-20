as yet unclassified
evaluateIteration
	| middleValue |
	result := positiveX * negativeXValue - (negativeX * positiveXValue) / (negativeXValue - positiveXValue).
	middleValue := functionBlock value: result.
	middleValue <= 0 ifTrue: [ negativeX := result. negativeXValue := middleValue. positiveXValue := 1 * positiveXValue ].
	middleValue >= 0 ifTrue: [ positiveX := result. positiveXValue := middleValue. negativeXValue := 1 * negativeXValue ].
	^ self relativePrecision: (positiveX - negativeX) abs