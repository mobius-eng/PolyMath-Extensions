instance creation
ewtRMSSolver: reference
	^ self new
		vectorNorm: (PMErrorWeightVectorNorm new
							defaultOfSize: reference size;
							reference: reference; yourself);
		desiredPrecision: 1.0;
		yourself.