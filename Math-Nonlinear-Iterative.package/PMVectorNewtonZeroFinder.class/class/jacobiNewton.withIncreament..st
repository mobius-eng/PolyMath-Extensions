instance creation
jacobiNewton: functionBlock withIncreament: increamentBlock
	| diagJacobian solver |
	solver := self new.
	diagJacobian := [ :y | | dy df |
		dy := increamentBlock value: y.
		df := (functionBlock value: y + dy) - solver functionValue.
		(df collectWithIndex: [ :dfn :n | dfn / (dy at: n) ]) asDhbVector ].
	solver setFunction: functionBlock; setJacobianOperator: diagJacobian; linearSolver: (PMDiagonalSolver new).
	^ solver