as yet unclassified
applyDirect: aVector
	^ ((1 to: aVector size) collect: [ :n |
		(diagonal at: n) * (aVector at: n) ]) asDhbVector