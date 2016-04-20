as yet unclassified
testLUPPlug
	| matrix solver residual |
	matrix := DhbMatrix new initializeRows: #(#(1 2 3) #(4 5 6) #(7 9 9)).
	solver := PMLUPPlug new linearOperator: matrix; rightHandSideVector: #(0 0 -2) asDhbVector.
	solver evaluate.
	self assert: solver hasConverged.
	residual := matrix * solver result - #(0 0 -2) asDhbVector.
	self assert: (PML1Norm new norm: residual) < DhbFloatingPointMachine new defaultNumericalPrecision.
	residual := solver result - #(1 -2 1) asDhbVector.
	self assert: (PML1Norm new norm: residual) < DhbFloatingPointMachine new defaultNumericalPrecision.