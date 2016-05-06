tests
testJacobiNewtonDiffusion
	"Solve discretized diffusion problem kD^2[c] - c^2 = 0
	with k = 0.1 on the interval (0,1) with mesh (0, 0.1, 0.2, ..., 1)
	and boundary conditions c(0)=1 and c(1) = 0.
	The model result is obtain via `fsolve` in Scilab"
	| solver modelAnswer |
	solver := PMVectorNewtonZeroFinder jacobiNewton: (self diffusionFunction: 0.1 dx: 0.1) withIncreament: [ :x | (x collect: [ :each | 0.1 * each + 1e-3 ]) asDhbVector ].
	solver initialValue: (1.0 to: 0.0 by: -0.1) asDhbVector.
	solver maximumIterations: 200.
	solver evaluate.
	modelAnswer := #(1.0 0.7790568 0.6188065 0.4968483 0.3995760 0.3182697 0.2470930 0.1820219 0.1202639 0.0599522 0.0) asDhbVector.
	self checkDiffusionSolution: solver with: modelAnswer.