private
setupDiffusionSolver: k dx: dx
	| solver |
	solver := PMVectorNewtonZeroFinder new.
	^ solver setFunction: (self diffusionFunction: k dx: dx);
		setJacobianOperator: (self diffusionJacobian: k dx: dx);
		initialValue: (1.0 to: 0 by: dx negated) asDhbVector;
		yourself