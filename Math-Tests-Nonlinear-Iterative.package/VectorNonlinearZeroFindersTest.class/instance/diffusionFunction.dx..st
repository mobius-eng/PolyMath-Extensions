private
diffusionFunction: k dx: dx
	^ [ :v | | y |
		y := DhbVector new: v size.
		y at: 1 put: (v at: 1) - 1.0; at: y size put: (v at: v size).
		2 to: y size - 1 do: [ :n |
			y at: n put: k * ((v at: n-1) - (2 * (v at: n)) + (v at: n+1)) /
								dx squared - (v at: n) squared ].
		y ]
