accessing
alpha: newAlpha
	"Sets alpha parameter of the algorithm:
	weight for derivative value in minimization criteria"
	(alpha isNumber and: [alpha > 0 & (alpha < 1)])
		ifFalse: [ self error: 'alpha is expected to be a number between 0 and 1. Given ', alpha asString ].
	alpha := newAlpha