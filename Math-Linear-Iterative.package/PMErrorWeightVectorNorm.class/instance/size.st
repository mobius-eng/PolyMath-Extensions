accessing
size
	relTolerance isNil ifFalse:[ ^relTolerance size ].
	absTolerance isNil ifFalse:[ ^absTolerance size ].
	reference isNil ifFalse:[ ^reference size. ].
	^0