initialization
setJacobianOperator: aBlock 
	"Defines the derivative of the function for which zeroes will be found."
	(aBlock respondsTo: #value:)
		ifFalse: [self error: 'Jacobian operator must implement the method value:value:'].
	"Trust that this block is correct derivative"
	jacobianOperator := aBlock