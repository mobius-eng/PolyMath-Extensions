operation
computeInitialValues
	initResidual := rightHandSideVector - (functionBlock value: result).
	residual := initResidual copy.
	residualChangeDirection := DhbVector new: residual size withAll: 0.
	searchDirection := residualChangeDirection copy.
	rho := 1.
	alpha := 1.
	omega := 1. 
	