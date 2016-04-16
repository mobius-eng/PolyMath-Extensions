initialization
setValueAtZero: valAtZero derivativeAtZero: derAtZero
	^ self
		setValueAtZero: valAtZero
		derivativeAtZero: derAtZero 
		valueAtOne: (functionBlock value: 1.0)
	