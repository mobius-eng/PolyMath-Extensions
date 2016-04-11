operation
residualError
	| tmp |
	tmp := 0.0.
	1 to: result size do: [ :n |
		tmp := ((residual at: n) / ((relTolerance at: n) * (reference at: n) abs + (absTolerance at: n))) squared + tmp. ].
	^ (tmp / result size) sqrt