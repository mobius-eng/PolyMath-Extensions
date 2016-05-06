*Math-Linear-Iterative
new: size withAll: aNumber
	"Creates square matrix with all elements equal to a number"
	| rows |
	rows := (1 to: size) collect: [ :n |
		DhbVector new: size withAll: aNumber ].
	^ self rows: rows