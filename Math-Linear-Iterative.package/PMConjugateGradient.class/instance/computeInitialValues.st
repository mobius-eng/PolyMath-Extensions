operation
computeInitialValues
	"r = b - A * x; p=r; beta = 0"
	residual := rightHandSideVector - (functionBlock value: result).
	residualSquared := residual * residual.
	searchDirection := residual copy.
	beta := 0.