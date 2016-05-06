instance creation
ewtRMSSolver: reference precision: precision
	^ self new
		vectorNorm: (PMErrorWeightVectorNorm new
							reference: reference;
							absTolerance: precision relTolerance: precision;
							yourself);
		desiredPrecision: 1.0;
		yourself.