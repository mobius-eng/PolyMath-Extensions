operation
rootMeanSquareNormOf: aVector withReference: ref
	| result |
	result := 0.
	absTolerance withIndexDo: [ :atol :i |
		result := ((aVector at: i) / ((relTolerance at: i) * (ref at: i) abs + atol)) squared + result ].
	result := result / (absTolerance size).
	^ result sqrt