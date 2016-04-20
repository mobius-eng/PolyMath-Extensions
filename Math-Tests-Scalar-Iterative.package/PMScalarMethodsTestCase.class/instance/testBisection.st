iterative algorithms
testBisection
	"Tests bisection on x*exp(x)-1=0 between 0 and 1"
	| solver |
	solver := PMBisectionZeroFinder function: [ :x | x * x exp - 1.0 ].
	solver setNegativeX: 0.0; setPositiveX: 1.0.
	solver evaluate.
	self assert: (solver result - 0.5671432909517471) abs <= DhbFloatingPointMachine new defaultNumericalPrecision.
	self assert: solver hasConverged = true.