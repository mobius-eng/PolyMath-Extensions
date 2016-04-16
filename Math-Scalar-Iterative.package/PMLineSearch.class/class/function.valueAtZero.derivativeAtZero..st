as yet unclassified
function: functionBlock valueAtZero: valAtZero derivativeAtZero: derAtZero
	^ self
		function: functionBlock
		valueAtZero: valAtZero
		derivativeAtZero: derAtZero
		valueAtOne: (functionBlock value: 1.0)