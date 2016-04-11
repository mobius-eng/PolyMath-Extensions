operation
computeInitialValues
	"Computes initial values as (1, g(1), 0, g(0))"
	result := 1.0.
	valueAtResult := valueAtOne.
	previousResult := 0.0.
	valueAtPreviousResult := valueAtZero.
	useCubicApproximation := false.
	(valueAtOne < valueAtZero or: [ valueAtOne <= minValue ])
			ifTrue: [ precision := 0.0. PMStopIterations new signal ].
