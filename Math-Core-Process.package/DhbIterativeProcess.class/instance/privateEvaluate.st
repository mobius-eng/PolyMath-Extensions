private
privateEvaluate
	"Evaluate internals: perform evaluations. Restart of stop iterations if signalled."
	[ 
		self initializeIterations.
		[
			iterations := iterations + 1.
			precision := self evaluateIteration.
			self hasConverged or: [ iterations >= maximumIterations ]
		] whileFalse
	] on: PMStopIterations do: [ :exception | ].
	self finalizeIterations.
	^self result