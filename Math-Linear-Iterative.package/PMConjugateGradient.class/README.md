I implement Conjugate Gradient method to solve the set of linear equations `A * x = b`. `A` can be either a matrix, or a linear operator:

```
A: vector -> vector
```

with properties `A(ax+by) = a A(x) + b A(y)`, where `x` and `y` are vectors  and `a` and `b` are scalars. `A` must be symmetric `A=A^t` and positively defined (`x^t * A * x > 0`)

## Example of use

```smalltalk
matrix := DhbMatrix rows: #(#(1 0 0) (0 2 0) (0 0 3)).
	rhsVector := #(1 2 3) asDhbVector.
	initValue := #(5 5 5) asDhbVector.
	cg := PMConjugateGradientL2 new.
	cg linearOperator: [ :x | matrix * x ]; rightHandSideVector: rhsVector; initialValue: initValue.
	cg desiredPrecision: DhbFloatingPointMachine new defaultNumericalPrecision.
	result := cg evaluate
```

Notice, `#linearOperator:` can take either a matrix (something that can be multiplied by a vector to produce a vector, i.e. responding to `#*`) or a function block (responding to `#value:`).