initialization
initialize
	super initialize.
	alpha := 1e-4.
	failingMin := 1e-3.
	"Since function is minimized along the line
	all we care for the precision to be negative"
	desiredPrecision := 0.0.
	minValue := DhbFloatingPointMachine new defaultNumericalPrecision.
	^ self