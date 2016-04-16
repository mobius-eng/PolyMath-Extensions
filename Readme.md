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

## Iterative linear solvers

[Conjugate gradient](https://en.wikipedia.org/wiki/Conjugate_gradient_method) method is an efficient method to solve sparse linear equations `A*x=b` with symmetric and positively defined matrix `A`.

`PMConjugateGradientL2` and `PMConjugateGradientEWT` provide implementations of unpreconditioned linear equations. The former estimates the residual error using [Euclidean (L2) norm](https://en.wikipedia.org/wiki/Euclidean_distance), whereas the latter, [root mean square norm with error weights](https://en.wikipedia.org/wiki/Root_mean_square).

`PMPreconditionedCGL2` implements (left) preconditioned CG method with Euclidean residual norm.

[Biconjugate gradient stabilised method  (BiCGStab)](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method) extends CG method for general matrices.

`PMBiCGStabL2` and `PMPreconditionedBiCGStabL2` implement BiCGStab with Euclidean norm, unpreconditioned and (left) preconditioned respectively.

## Future plans

These extensions form a part of my PhD project, thus I'm implementing what is necessary for it. The ultimate goal is to implement efficient ODE solver in Smalltalk, similar to LSODA or ODEPACK. Iterative linear solvers are the first steps to this way. The following steps will be taken further:

- Nonlinear solvers: Newton's method for set of equations, Broyden's method.
- Predictor-corrector ODE steppers (with iterative corrector) based on Adams-Moulton (AM) and backward-differences formulas (BDF).
- Error controlled steppers with step and method order controls.
- Automatic switch between AM and BDF methods.

## Clashes with PolyMath packages

To minimize the clashes with PolyMath and its updates, almost all the classes are located in packages that are not present in PolyMath. Exception is made for `DhbIterativeProcess` since it is the original class that is affected.
