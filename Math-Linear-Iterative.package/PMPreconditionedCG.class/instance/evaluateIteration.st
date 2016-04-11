operation
evaluateIteration
	| alpha newResidualSquared |
	"p = z + beta * p"
	1 to: result size do: [ :n | searchDirection at: n put: (beta * (searchDirection at: n) + (conditionedResidual at: n)) ].
	"ap = A * p"
	residualChangeDirection := functionBlock value: searchDirection.
	"alpha = (z,r) / (p, ap)"
	alpha := residualSquared  / (searchDirection * residualChangeDirection).
	"In one go:
		x = x + alpha * p
		r = r - alpha * ap"
	1 to: result size do: [ :n |
		result at: n put: (result at: n) + (alpha * (searchDirection at: n)).
		residual at: n put: (residual at: n) - (alpha * (residualChangeDirection at: n)). ].
	conditionedResidual := preconditioner applyInverse: residual.
	newResidualSquared := conditionedResidual * residual.
	"beta = (r-new, r-new) / (r,r)"
	beta := newResidualSquared / residualSquared.
	residualSquared := newResidualSquared.
	^ self residualError