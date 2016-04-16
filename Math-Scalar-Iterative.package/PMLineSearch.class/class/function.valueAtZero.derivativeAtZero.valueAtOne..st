as yet unclassified
function: functionBlock valueAtZero: valAtZero derivativeAtZero: derAtZero valueAtOne: valAtOne
	| instance |
	instance := self new.
	instance setFunction: functionBlock;
		setValueAtZero: valAtZero derivativeAtZero: derAtZero valueAtOne: valAtOne.
	^ instance