private
vectorize: anObject
	(anObject isNumber and: [ self size > 0 ])
		ifTrue: [ ^DhbVector new: (self size) withAll: anObject ].
	(anObject isCollection)
		ifTrue: [ 
			(self size = 0 or: [ self size = anObject size ])
				ifTrue: [ ^anObject asDhbVector ]
				ifFalse: [ self error: 'Size mismatch' ] ].
	self error: 'Unknown object'