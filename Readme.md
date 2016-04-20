# PolyMath - extensions
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/PolyMathOrg/PolyMath/master/LICENSE)

[*PolyMath*](https://github.com/PolyMathOrg/PolyMath) is mathematical library for Smalltalk ([Pharo](http://pharo.org/)). As the library is in the transition state, this repository consists of some experimental features that are in the process of being integrated with PolyMath. **Highly experimental: use on your own risk.**

## Iterative process

*PolyMath* implements iterative algorithms using generic iterative process construct encapsulated in class `DbhIterativeProcess`. While this class is powerful to represent many iterative processes (Newton's method, functional iterations), it is not straight-forward to end the process before it has converged *and* before the maximum iterations was reached. For example, it might became apparent that the process wouldn't converge in the middle of iteration. Or, in slightly different case, the convergence might need to be checked in the middle of the single iteration and if it is reached, it would be desirable to exit unconditionally (and spare one more convergence check at the end).

The implementation provided by this repository introduces exception `PMStopIterations`, which can be raised by a particular implementation. Once signaled, the iterative process will be finished. This extension is back-compatible with all dependent algorithms in *PolyMath*. This solution might look like a heavy hack, as the use of exceptions for anything but errors is not advised. So far, considering multiple constraints, it is the best solution found, but the search for a cleaner approach continues.

In addition, it is deliberated now if an unhandled exception *should* be raised if the process fails to converge. This, however, would almost certainly break back-compatibility.

## Iterative scalar solvers
*PolyMath* already contains line search method (under `DhbLineSearch`, authored by mobius-eng) that finds the minimum of the function between 0 and 1. However, it does contain a couple of small bugs that are really hard to pick up in tests. But it is difficult to reliably update it in *PolyMath* directly due to transition and massive updates that it undergoes. So, for now, this repository offers package `Math-Scalar-Iterative` that contains updated version of this algorithm.

In addition, it also contains globally convergent scalar Newton's method `PMScalarNewtonZeroFinder` (it uses line search to avoid usual pitfalls of unmodified Newton's method). The name deliberately contains "scalar" to distinguish it from more general vector-based Newton's method for the set of equations.

Other solvers: `PMBisectionZeroFinder` and `PMFalsePositionZeroFinder`. The former is almost identical to `DhbBisectionZeroFinder` except it treats properly th case when the exact solution was found during the iteration (extremely rare, but possible case, for example solving `x squared - 4 = 0` on `(1, 3)`). The latter implements modified *false position* method (see Hamming R.W., *Numerical methods for scientists and engineers* for more details).

## Iterative linear solvers

**There was a substantial refactoring** performed on linear iterative solvers. The residual norm evaluation was delegated to an instance variable `vectorNorm` (must respond to `#norm:` message with the vector as an argument). Four possible norms were introduced: `PML1Norm`, `PML2Norm`, `PMMaxNorm` and `PMErrorWeightVectorNorm` that implement l1-, l2-, l-infinity- and root-mean-square (RMS) with error weights vector (EWT) norms. L1-norm is a sum of absolute values of the vector. [l2 (Euclidean) norm](https://en.wikipedia.org/wiki/Euclidean_distance) takes square root of sum of squares. L-infinity takes the maximum of absolute values. [Root mean square norm with error weights](https://en.wikipedia.org/wiki/Root_mean_square) uses more complex formula in case vector components are heterogeneous.

[Conjugate gradient](https://en.wikipedia.org/wiki/Conjugate_gradient_method) method is an efficient method to solve sparse linear equations `A*x=b` with symmetric and positively defined matrix `A`.

`PMConjugateGradient` and `PMPreconditionedCG` provide the implementations for unpreconditioned and (left) preconditioned linear equations respectively.

[Biconjugate gradient stabilised method  (BiCGStab)](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method) extends CG method for general matrices.

`PMBiCGStab` and `PMPreconditionedBiCGStab` implement BiCGStab unpreconditioned and (left) preconditioned respectively.

On the top of this, to provide the uniformity between the solvers, `PMLUPPlug` is an adapter of an existing LUP-decomposition solver for matrices `DhbLUPDecomposition`.

## Nonlinear solvers

Globally convergent vector-based Newton's method is implemented by `PMVectorNewtonZeroFinder` in package `Math-Nonlinear-Iterative`.

Couple of tests for this method are provided, but it has not been rigorously tested yet! **Use at your own risk**.

## Future plans

These extensions form a part of my PhD project, thus I'm implementing what is necessary for it. The ultimate goal is to implement efficient ODE solver in Smalltalk, similar to LSODA or ODEPACK. Iterative linear solvers are the first steps to this way. The following steps will be taken further:

- Nonlinear solvers: Broyden's method.
- Predictor-corrector ODE steppers (with iterative corrector) based on Adams-Moulton (AM) and backward-differences formulas (BDF).
- Error controlled steppers with step and method order controls.
- Automatic switch between AM and BDF methods.

On the top of the priority list is to provide a proper documentation and examples for the methods. It will be likely provided as a Github project site (just because it can deal with mathematical formulas better than plain Markdown).

## Clashes with PolyMath packages

To minimize the clashes with PolyMath and its updates, almost all the classes are located in packages that are not present in PolyMath. Exception is made for `DhbIterativeProcess` since it is the original class that is affected.
