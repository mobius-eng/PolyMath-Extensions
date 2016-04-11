I am the root for Linear Iterative methods

Implemented methods:

`initialValue:` sets initial approximation
`linearOperator:` sets  linear operator (in form of a function or a matrix)
`residual` answers residual vector
`rightHandSideVector:` sets the right hand side vector `b` for `Ax=b`.

Abstract methods

`residualError` answers the error (number) of the residual, i.e. how close the result is to actual solution.

There are multiple ways to calculate it: using Euclidean (L2) norm, maximum abs element (L1 norm), or more
elaborate, root-mean-square error using error-weight vector. Subclasses must implement a particular way.

Inhereted abstract methods

Subclasses still must implement

`computeInitialValues` and `evaluateIteration`
