---
layout: post
title: "Line Search Algorithm"
date: "2016-05-06 14:45:13 +0200"
categories: numerical-methods scalar
---

## Line Search Algorithm

It is often happens that minimization algorithms (including zero finders, such as Newton's method or Broyden's method) fail to converge if initial approximation was chosen far away from the actual solution. Nevertheless, it is possible to improve the current approximation in the direction of minimization even if the method suggests too large a step that overshoots the solution. Consider trivial example of solving $$\operatorname{atan}{x}=0$$ with solution $$x=0$$. Newton's algorithm will result in the following recursive formula:

$$
x_{n+1} = x_n - \frac{\operatorname{atan}{x_n}}{1/(1+x_n^2)}
$$

If the initial approximation $$ | x_0 | \geq 1.4$$, for example $$x_0 = 2$$, the sequence will diverge:

$$
x_0 = 2, x_1 = -3.54, x_2 = 13.95, x_3 = -279.3,\ldots
$$

However, it possible to recover from this divergence since Newton's method takes the step towards the minimization of $$f(x) =  \operatorname{atan}^2{x} $$, but overshoots. Thus, it is reasonable to find a better approximation somewhere between $$x_0$$ and $$x_1$$. One way to do this is to use bisection method. Indeed, we know the minimum exists between these two points. So, the first point to try would be $$x'=(x_0+x_1)/2$$. A different approach was suggested in *Numerical Recipes* by W. H. Press *et al*: we usually know $$f(x_0)$$, $$D[f](x_0)$$ (derivative of $$f(x)$$ at $$x_0$$) and $$f(x_1)$$ --- these points can be used to approximate $$f(x)$$ on the line between $$x_0$$ and $$x_1$$ via quadratic function.
