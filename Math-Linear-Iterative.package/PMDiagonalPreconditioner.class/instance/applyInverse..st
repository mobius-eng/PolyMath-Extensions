as yet unclassified
applyInverse: aVector
	^ ((1 to: aVector size) collect: [ :n |
		(inverseDiagonal at: n) * (aVector at: n) ]) asDhbVector