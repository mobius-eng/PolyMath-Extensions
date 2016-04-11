as yet unclassified
evaluateIteration
	| t beta newRho |
	newRho := initResidual * residual.
	"Restart if newRho became too small"
	newRho abs < (DhbFloatingPointMachine new defaultNumericalPrecision * rhsSquared)
		ifTrue: [ PMRestartIterations new signal ].
	beta := (newRho / rho) * (alpha / omega).
	rho := newRho.
	searchDirection := beta * (searchDirection - (omega * residualChangeDirection )) + residual.
	residualChangeDirection := functionBlock value: searchDirection.
	alpha := rho / (initResidual * residualChangeDirection).
	result withIndexDo: [ :x :n |
		result at: n put: (x + (alpha * (searchDirection at: n))).
		residual at: n put: ((residual at: n) - (alpha * (residualChangeDirection at: n))) ].
	"Finish if residual got small"
	precision := self residualError.
	self hasConverged ifTrue: [ ^ precision ].
	"Extra bit if residual is still large"
	t := functionBlock value: residual.
	omega := (t * residual) / (t * t).
	result withIndexDo: [ :x :n |
		result at: n put: (x + (omega * (residual at: n))).
		residual at: n put: ((residual at: n) - (omega * (t at: n))) ].
	^ self residualError