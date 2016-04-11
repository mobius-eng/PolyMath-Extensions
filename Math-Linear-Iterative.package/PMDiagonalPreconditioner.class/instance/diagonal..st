as yet unclassified
diagonal: aVector
	diagonal := aVector asDhbVector.
	inverseDiagonal := (diagonal collect: [ :each | each reciprocal ]) asDhbVector