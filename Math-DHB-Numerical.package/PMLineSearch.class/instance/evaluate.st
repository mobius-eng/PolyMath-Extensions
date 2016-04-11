operation
evaluate
	"Overrides the parent method to check if method needs to be executed or
	the full step (x=1) does minimize the function already"
	(valueAtOne < valueAtZero
		or: [ valueAtOne <= minValue ])
			ifTrue: [^ result := 1. ].
	^ super evaluate