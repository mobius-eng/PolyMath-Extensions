private
checkDiffusionSolution: solver with: modelResult
	| diff |
	self assert: solver hasConverged.
	diff := modelResult asDhbVector - solver result.
	"The answer is known within 1e-7, thus have to relax the tolerance"
	self assert: solver hasConverged.
	self assert: (PMErrorWeightVectorNorm new defaultOfSize: 11; absTolerance:1.0e-7 relTolerance:1.0e-7; norm: diff ) <= 1.0.
	self assert: (PMErrorWeightVectorNorm new defaultOfSize: 11; norm: (solver functionValue)) <= 1.0