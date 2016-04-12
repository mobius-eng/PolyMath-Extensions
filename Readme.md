# PolyMath - extensions

[*PolyMath*](https://github.com/PolyMathOrg/PolyMath) is mathematical library for Smalltalk ([Pharo](http://pharo.org/)). As the library is in the transition state, this repository consists of some experimental features that are in the process of being integrated with PolyMath. **Highly experimental: use on your own risk.**

## Iterative process

*PolyMath* implements iterative algorithms using generic iterative process construct encapsulated in class `DbhIterativeProcess`. While this class is powerful to represent many iterative processes (Newton's method, functional iterations), it is not straight-forward to end the process before it has converged *and* before the maximum iterations was reached (for example, when it became apparent that the process wouldn't converge).

The implementation provided by this repository introduces exception `PMStopIterations`, which can be raised by a particular implementation. Once signaled, the iterative process will be finished. This extension is back-compatible with all dependent algorithms in *PolyMath*.

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
