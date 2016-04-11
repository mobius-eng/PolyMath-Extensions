accessing
failingMin: newFailingMin
	"Sets abs minimum of the step, below which the method will produce failure"
	(newFailingMin isNumber and: [ (newFailingMin < 1) & (newFailingMin >= 0) ])
		ifFalse: [
			self error: 'failingMin must be a number between 0 and 1. Given ', newFailingMin asString ].
	failingMin := newFailingMin