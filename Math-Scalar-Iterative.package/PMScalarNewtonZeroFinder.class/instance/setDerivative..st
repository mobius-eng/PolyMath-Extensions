initialization
setDerivative: aBlock 
	"Defines the derivative of the function for which zeroes will be found."
	(aBlock respondsTo: #value:)
		ifFalse: [self error: 'Derivative block must implement the method value:'].
	"Trust that this block is correct derivative"
	derivativeBlock := aBlock