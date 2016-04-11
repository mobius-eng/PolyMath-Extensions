operation
computeInitialValues
	"Computes initial values as (1, g(1), 0, g(0))"
	result := 1.0.
	valueAtResult := valueAtOne.
	previousResult := 0.0.
	valueAtPreviousResult := valueAtZero.
	useCubicApproximation := false.
