private
diffusionJacobian: k dx: dx
	^ [ :c |
		| m |
		m := DhbMatrix new: c size withAll: 0.0.
		m at: 1 at: 1 put: 1.0.
		m at: c size at: c size put: 1.0.
		2 to: c size - 1 do: [ :n |
			m at:n at:n put: -2.0 * k / dx squared - (2.0 * (c at: n)).
			m at:n at:n-1 put: k / dx squared.
			m at:n at:n+1 put: k / dx squared ].
		m ]