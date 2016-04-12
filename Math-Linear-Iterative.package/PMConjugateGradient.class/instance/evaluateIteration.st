operation
evaluateIteration
	| alpha newResidualSquared |
	"p = r + beta * p"
	1 to: result size do: [ :n |
		searchDirection at: n put: (beta * (searchDirection at: n) + (residual at: n)) ].
	"ap = A * p"
	residualChangeDirection := functionBlock value: searchDirection.
	"alpha = (r,r) / (p, ap)"
	alpha := residualSquared  / (searchDirection * residualChangeDirection).
	"In one go:
		x = x + alpha * p
		r = r - alpha * ap
		and new_rr = (r,r)"
	newResidualSquared := 0.0.
	1 to: result size do: [ :n |
		result at: n put: (result at: n) + (alpha * (searchDirection at: n)).
		residual at: n put: (residual at: n) - (alpha * (residualChangeDirection at: n)).
		newResidualSquared := newResidualSquared + (residual at: n) squared. ].
	"beta = (r-new, r-new) / (r,r)"
	beta := newResidualSquared / residualSquared.
	residualSquared := newResidualSquared.
	^ self residualError