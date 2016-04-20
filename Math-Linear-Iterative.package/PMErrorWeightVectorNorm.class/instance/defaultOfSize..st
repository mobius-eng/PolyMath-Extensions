as yet unclassified
defaultOfSize: anInteger
	^ self uniformSize: anInteger
		absTolerance: DhbFloatingPointMachine new defaultNumericalPrecision
		relTolerance: DhbFloatingPointMachine new defaultNumericalPrecision
	 