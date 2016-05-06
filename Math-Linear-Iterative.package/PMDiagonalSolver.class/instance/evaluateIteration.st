as yet unclassified
evaluateIteration
	result := (rightHandSideVector collectWithIndex: [ :b :n |
		b / (functionBlock at: n) ]) asDhbVector.
	^ 0.0