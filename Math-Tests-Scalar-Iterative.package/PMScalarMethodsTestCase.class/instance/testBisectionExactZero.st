iterative algorithms
testBisectionExactZero
	"Tests successful return from bisection if exact zero is found"
	| solver |
	solver := PMBisectionZeroFinder function: [ :x | x squared - 4 ].
	solver setNegativeX: 1.0; setPositiveX: 3.0.
	solver evaluate.
	self assert: solver result = 2.0.
	self assert: solver iterations = 1.
	self assert: solver hasConverged = true.