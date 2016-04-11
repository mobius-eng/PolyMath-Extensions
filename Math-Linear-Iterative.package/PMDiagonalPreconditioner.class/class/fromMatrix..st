as yet unclassified
fromMatrix: aMatrix
	| instance  dimension diagonal |
	dimension := aMatrix dimension.
	diagonal := ((1 to: dimension x) collect: [ :i | aMatrix at: i at: i ]) asDhbVector.
	instance := self new.
	instance diagonal: diagonal.
	^ instance