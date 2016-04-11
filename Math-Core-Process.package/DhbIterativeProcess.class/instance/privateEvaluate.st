private
privateEvaluate
	"Evaluate internals: perform evaluations. Restart of stop iterations if signalled."
	[ 
		self initializeIterations.
		[
			iterations := iterations + 1.
			precision := self evaluateIteration.
			self hasConverged or: [ iterations >= maximumIterations ]
		] whileFalse: []
	] on: PMRestartIterations,PMStopIterations do: [ :exception |
		(exception isKindOf: PMStopIterations) ifTrue: [  ].
		(exception isKindOf: PMRestartIterations) ifTrue: [
			iterations < maximumIterations ifTrue: [ ^ self privateEvaluate ] ] ].
	self finalizeIterations.
	^self result