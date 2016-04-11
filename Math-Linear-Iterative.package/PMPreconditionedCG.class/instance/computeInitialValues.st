operation
computeInitialValues
	super computeInitialValues.
	conditionedResidual := preconditioner applyInverse: residual.
	searchDirection := conditionedResidual copy.
	residualSquared := conditionedResidual * residual.