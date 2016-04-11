operation
evaluateIteration
	| a b tmp1 tmp2 gamma1 gamma2 nextX deltaX |
	deltaX := result - previousResult.
	useCubicApproximation
		ifFalse: [
			nextX := derivativeAtZero negated * 0.5 / (valueAtOne - valueAtZero - derivativeAtZero) max: 0.1.
			useCubicApproximation := true ]
		ifTrue: [
			gamma1 := valueAtResult - (derivativeAtZero * result) - valueAtZero / deltaX.
			gamma2 := valueAtPreviousResult - (derivativeAtZero * previousResult) - valueAtZero / deltaX.
			tmp1 := gamma1 / result squared.
			tmp2 := gamma1 / previousResult squared.
			a := tmp1 - tmp2.
			b := result * tmp2 - previousResult * tmp1.
			nextX := (b negated + (b squared - (3.0 * a * derivativeAtZero) sqrt)) / (3.0 * a) min: 0.5 * result max: 0.1 * result.
			nextX < failingMin ifTrue: [ self error: 'next x is too small: cannot minimize' ] ].
	self updateResult: nextX.
	^ valueAtResult - (alpha * derivativeAtZero + valueAtZero)