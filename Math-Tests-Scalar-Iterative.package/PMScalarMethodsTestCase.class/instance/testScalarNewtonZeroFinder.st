iterative algorithms
testScalarNewtonZeroFinder
	"Code Example 5.3"

	| zeroFinder result |
	zeroFinder := PMScalarNewtonZeroFinder
		function: [ :x | x errorFunction - 0.9 ]
		derivative: [:x | 2.0 / Float pi sqrt * x squared negated exp ].
	zeroFinder initialValue: 1.0.
	result := zeroFinder evaluate.
	self assert: zeroFinder hasConverged.
	self assert: (result - 1.28155193867885) abs < zeroFinder precision