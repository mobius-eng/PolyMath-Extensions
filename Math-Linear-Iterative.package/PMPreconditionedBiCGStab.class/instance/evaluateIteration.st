as yet unclassified
evaluateIteration
	| y z t u beta newRho |
	newRho := initResidual * residual.
	"Restart if newRho became too small"
	newRho abs < (DhbFloatingPointMachine new defaultNumericalPrecision * rhsSquared)
		ifTrue: [ self initializeIterations. ^ self evaluateIteration ].
	beta := (newRho / rho) * (alpha / omega).
	rho := newRho.
	searchDirection := beta * (searchDirection - (omega * residualChangeDirection )) + residual.
	y := preconditioner applyInverse: searchDirection.
	residualChangeDirection := functionBlock value: y.
	alpha := rho / (initResidual * residualChangeDirection).
	result withIndexDo: [ :x :n |
		result at: n put: (x + (alpha * (y at: n))).
		residual at: n put: ((residual at: n) - (alpha * (residualChangeDirection at: n))) ].
	"Finish if residual got small"
	precision := self residualError.
	self hasConverged ifTrue: [ ^ precision ].
	z := preconditioner applyInverse: residual.
	"Extra bit if residual is still large"
	t := functionBlock value: z.
	u := preconditioner applyInverse: t.
	omega := (u * z) / (u * u).
	result withIndexDo: [ :x :n |
		result at: n put: (x + (omega * (z at: n))).
		residual at: n put: ((residual at: n) - (omega * (t at: n))) ].
	^ self residualError