initialization
uniformSize: anInteger absTolerance: atol relTolerance: rtol
	absTolerance := DhbVector new: anInteger withAll: atol.
	relTolerance := DhbVector new: anInteger withAll: rtol.
	reference := DhbVector new: anInteger withAll: 1.0.
	^ self