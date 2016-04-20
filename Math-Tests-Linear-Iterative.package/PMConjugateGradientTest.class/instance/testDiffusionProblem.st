tests
testDiffusionProblem
	"Tests a diffusion-like problem solution. The linear operator is defined as a function"
	| diffusionOperator rhsVector initValue cg result solution |
	"slightly faux: the first and last elements of the result are meaningless"
	diffusionOperator := [ : v |
		| diff deltaX w |
		diff := 1.0.
		deltaX := 1.0.
		w := DhbVector new: v size.
		w at: 1 put: ((v at: 1) + (v at:2)).
		2 to: v size-1 do: [ :n |
			w
				at: n
				put: (diff * ((v at: (n-1)) - (2 * (v at: n)) + (v at: (n+1))) / deltaX squared) ].
		w at: w size put: ((v at: v size) + (v at: v size - 1)).
		w ].
	rhsVector := #(3 0 0 6) asDhbVector.
	initValue := #(1 1 1 1) asDhbVector.
	cg := PMConjugateGradient ewtRMSSolverFor: diffusionOperator rightHandSide: rhsVector initialValue: initValue.
	result := cg evaluate.
	solution := (diffusionOperator value: result) - rhsVector.
	self assert: cg hasConverged.
	self assert: (PMErrorWeightVectorNorm new defaultOfSize: rhsVector size; reference: rhsVector; norm: solution) <= 1.