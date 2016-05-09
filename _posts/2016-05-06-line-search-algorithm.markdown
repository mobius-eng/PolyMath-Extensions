---
layout: post
title: "Line Search Algorithm"
date: "2016-05-06 14:45:13 +0200"
categories: numerical-methods scalar
---

## Line Search Algorithm

It is often happens that minimization algorithms (including zero finders, such as Newton's method or Broyden's method) fail to converge if initial approximation was chosen far away from the actual solution. Nevertheless, it is possible to improve the current approximation in the direction of minimization even if the method suggests too large a step that overshoots the solution. Consider trivial example of solving $$\operatorname{atan}{x}=0$$ with solution $$x=0$$ using Newton's method:

$$
x_{n+1} = x_n - \frac{\operatorname{atan}{x_n}}{1/(1+x_n^2)}
$$

If the initial approximation $$ \vert x_0\vert\geq1.4$$, for example $$x_0 = 2$$, the sequence will diverge:

$$
x_0 = 2, x_1 = -3.54, x_2 = 13.95, x_3 = -279.3,\ldots
$$

However, it possible to recover from this divergence since Newton's method takes the step towards the minimization of $$f(x) =\vert\operatorname{atan}{x}\vert$$, but overshoots. Thus, it is reasonable to find a better approximation somewhere between $$x_0$$ and $$x_1$$. One way to do this is to use bisection method. Indeed, we know the minimum exists between these two points. So, the first point to try would be $$x'=(x_0+x_1)/2$$ (see R.W. Hamming *Numerical methods for scientists and engineers*). A different approach was suggested in *Numerical Recipes* by W.H. Press *et al*: we usually know $$f(x_0)$$, $$D[f](x_0)$$ (derivative of $$f(x)$$ at $$x_0$$) and $$f(x_1)$$ --- these points can be used to approximate $$f(x)$$ on the line between $$x_0$$ and $$x_1$$ via quadratic function. Since this method is also applicable for vector-valued problem $$\boldsymbol{F}(\boldsymbol{x})=\boldsymbol{0}$$ for $$\boldsymbol{F}:\mathbb{R}^n\to\mathbb{R}^n$$, the following derivation will be conducted for a more general vector-valued problem.


### Derivation

Consider the iteration of a non-linear solver that suggests to take the step $$\boldsymbol{p} = \boldsymbol{x}_1 - \boldsymbol{x}_0$$. The line search algorithm aims to find such $$\lambda$$ that $$\boldsymbol{x}_0+\lambda\boldsymbol{p}$$ will minimize the absolute value of $$\boldsymbol{F}(\boldsymbol{x})$$. This can be stated as a minimization of $$f(\boldsymbol{x})=1/2\boldsymbol{F}(\boldsymbol{x})\cdot\boldsymbol{F}(\boldsymbol{x})$$, $$f:\mathbb{R}^n\to\mathbb{R}$$. Since it aims to find the minimum on the line $$(\boldsymbol{x}_0,\boldsymbol{x}_1)$$, we can simplify this problem further by introducing $$g(\lambda)=f(\boldsymbol{x}_0+\lambda\boldsymbol{p})$$ and minimizing it for $$\lambda\in(0,1)$$. Derivative of $$g(\lambda)$$:

$$
D[g](\lambda) = \boldsymbol{F}(\boldsymbol{x}_0+\lambda\boldsymbol{p})\cdot \mathbf{J}(\boldsymbol{x}+\lambda\boldsymbol{p})\cdot\boldsymbol{p}
$$

where $$\mathbf{J}(\boldsymbol{x})$$ is the Jacobian matrix $$\vert\vert \partial_j F_i \vert\vert$$. When $$\lambda=0$$ this simplifies to $$D[g](0) = \boldsymbol{F}(\boldsymbol{x}_0)\cdot \mathbf{J}(\boldsymbol{x})\cdot\boldsymbol{p}$$. Notice, if $$\boldsymbol{p}$$ is a Newtonian step: $$\boldsymbol{p} = - \mathbf{J}^{-1}(\boldsymbol{x}_0)\cdot \boldsymbol{F}(\boldsymbol{x}_0)$$, it simplifies further to $$D[g](0) = -2g(0)$$.

Before introducing iteration steps it should be noted that the line search only needs to be conducted if $$g(1)\leq g(0)$$. In case the line search is required, $$g(1) > g(0)$$, the stop criteria is chosen as:

$$
g(\lambda)\leq g(0) + \alpha \lambda D[g](0)
$$

where $$\alpha$$ is an empirical parameter (defaults to $$10^{-4}$$).

To find the minimum of $$g(x)$$, it can be approximated as a quadratic function using $$g(0)$$, $$D[g](0)$$ and $$g(1)$$ with minimum found at

$$
\lambda^* =- \frac{D[g](0)}{2[g(1)-g(0)-D[g](0)]}
$$

Since $$g(1)>g(0)$$, $$\lambda^* \lessapprox 1/2$$. If $$\lambda^* $$ satisfies minimization criteria, the method stops. Otherwise, $$g(x)$$ can be further approximated using cubic approximation using $$g(0)$$, $$D[g](0)$$ and $$g(\lambda)$$ values on last two approximations starting with $$\lambda_2=1$$ (second last) and $$\lambda_1=\lambda^* $$ (last):

$$
g(x)=ax^3+bx^2+D[g](0)x+g(0)
$$

Coefficients $$a$$ and $$b$$ found as:

$$
\gamma(\lambda) = \frac{1}{\lambda_1 - \lambda_2}\frac{1}{\lambda}[g(\lambda) - D[g](0)\lambda - g(0)]
$$

$$
a = \gamma(\lambda_1) - \gamma(\lambda_2)
$$

$$
b = -\lambda_2\gamma(\lambda_1) + \lambda_1\gamma(\lambda_2)
$$

And minimum of the cubic function is at:

$$
\lambda = \frac{-b + \sqrt{b^2 - 3aD[g](0)}}{3a}
$$

Finally, to prevent $$\lambda$$ to be too large or too small, at each step it is forced to fall in $$(0.1\lambda_1, 0.5\lambda_1)$$. If $$\lambda$$ becomes too small for the method to converge (say, below $$10^{-5}$$), the method fails.

### Code (with annotations)

In PolyMath-Extensions line search is implemented as a class `PMLineSearchFunctionMinimizer` in package `Math-Scalar-Iterative`. This class is descendant of `DhbIterativeProcess` (`PMIterativeProcess`), thus, main functionality is implemented in the method `evaluateIteration`:

```smalltalk
PMLineSearchFunctionMinimizer>>evaluateIteration
	| a b tmp1 tmp2 gamma1 gamma2 nextX deltaX |
	deltaX := result - previousResult. "λ1 - λ2"
	useCubicApproximation
		ifFalse: [
			"First step - quadratic approximation"
			"λ* = D[g](0)/2(g(1) - g(0) - D[g](0))"
			nextX := derivativeAtZero negated * 0.5 / (valueAtOne - valueAtZero - derivativeAtZero) max: 0.1.
			"After first step - use cubic approximation"
			useCubicApproximation := true ]
		ifTrue: [
			gamma1 := valueAtResult - (derivativeAtZero * result) - valueAtZero / deltaX.
			gamma2 := valueAtPreviousResult - (derivativeAtZero * previousResult) - valueAtZero / deltaX.
			tmp1 := gamma1 / result squared.
			tmp2 := gamma2 / previousResult squared.
			"a and b coefficients"
			a := tmp1 - tmp2.
			b := result * tmp2 - (previousResult * tmp1).
			nextX := (b negated + (b squared - (3.0 * a * derivativeAtZero) sqrt)) / (3.0 * a) min: 0.5 * result max: 0.1 * result.
			"Finish (fail) iterations if next approximation is too small"
			nextX < failingMin ifTrue: [ PMStopIterations new signal ] ].
	"Update λ2 ← λ1, λ1 ← nextX and update function values"
	self updateResult: nextX.
	"precision in this case: just need to be negative"
	^ valueAtResult - (alpha * derivativeAtZero + valueAtZero)
```

This class also overrides `computeInitialValues` and stops the evaluation if $$g(1)<g(0)$$:

```smalltalk
PMLineSearchFunctionMinimizer>>computeInitialValues
	"Computes initial values as (1, g(1), 0, g(0))"
	result := 1.0.
	valueAtResult := valueAtOne.
	previousResult := 0.0.
	valueAtPreviousResult := valueAtZero.
	useCubicApproximation := false.
	(valueAtOne < valueAtZero or: [ valueAtOne <= minValue ])
			ifTrue: [ precision := 0.0. PMStopIterations new signal ].
```

### Example of use

Line search is not intended to be used directly by the user, but rather by an implementer of other methods. Direct use examples can be found in test-case class `PMScalarMethodsTestCase`. Here is the example of its use in vector-vales Newton's method:

```smalltalk
PMVectorNewtonZeroFinder>>evaluateIteration
	"Compute one step of Newton's zero finding method. Answers the estimated precision."
	| g0 dg0 coefficient |
	"Calculate new Newton step"
	lastFunctionValue := newFunctionValue.
	linearSolver linearOperator: (jacobianOperator value: result);
		rightHandSideVector: lastFunctionValue; initialValue: searchStep negated.
	searchStep := linearSolver evaluate negated.
	linearSolver hasConverged ifFalse: [ PMStopIterations new signal ].
	"Initialize variables for line search"
	g0 := 0.5 * (lastFunctionValue * lastFunctionValue).
	dg0 := -2 * g0 squared.
	"Get coefficient of the step from line search: 0 < coefficient <= 1; newFunctionValue is updated as well"
	coefficient := lineSearch setValueAtZero: g0 derivativeAtZero: dg0; evaluate.
	lineSearch hasConverged ifFalse: [ PMStopIterations new signal ].
	result := coefficient * searchStep + result.
	^ searchStep rootMeanSquareNormWith: errorWeightVector
```

Quite importantly, this implementation uses the fact that line search guarantees that the nonlinear function $$\boldsymbol{F}(\boldsymbol{x})$$ will be called on $$\boldsymbol{x}_\text{new}=\boldsymbol{x}+\lambda\boldsymbol{p}$$, thus, avoiding unnecessary function evaluations.
