# Weierstrass Approximation via Bernstein Polynomials

The question is: can every continuous real-valued function on $[0,1]$ be approximated arbitrarily well by a polynomial? In mathematical language, given a continuous function $f:[0,1]\to\mathbb{R}$ and any tolerance $\varepsilon>0$, can we find a polynomial $p$ for which

\[
\sup_{x\in[0,1]} |f(x)-p(x)|<\varepsilon?
\]

Why care? Polynomials are unusually simple objects. A four-function calculator can evaluate them, their derivatives and integrals are immediate, and many algorithms are easier to analyze or optimize when the input is polynomial.

But hasn't calculus already solved this with Taylor polynomials? **Wrong.** Taylor approximation depends on derivative information at a point. The Weierstrass theorem assumes only continuity not differentiability.

In fact, in the precise sense of Baire category, a typical continuous real-valued function is nowhere differentiable. Such functions behave more like the Weierstrass function below than like the smooth curves used in introductory calculus. The explorer plots partial sums of the classical example

\[
W(x)=\sum_{k=0}^{\infty}0.6^k\cos(11^k\pi x).
\]

<iframe
  src="./assets/weierstrass-function.html"
  title="Interactive explorer for a Weierstrass function"
  width="100%"
  height="520"
  loading="lazy"
  sandbox="allow-scripts"
></iframe>

[Open the Weierstrass function explorer](./assets/weierstrass-function.html)

Taylor's method is still instructive. It anchors the approximation at one point and uses more derivative information as the degree increases. For cosine, centered at zero, the approximations are

\[
T_{2m}(x)=\sum_{j=0}^{m}\frac{(-1)^j x^{2j}}{(2j)!}.
\]

<iframe
  src="./assets/cosine-taylor-convergence.html"
  title="Taylor polynomials converging to cosine"
  width="100%"
  height="520"
  loading="lazy"
  sandbox="allow-scripts"
></iframe>

[Open the cosine Taylor explorer](./assets/cosine-taylor-convergence.html)

Bernstein's construction needs no derivatives at all. It samples the function on a grid and uses probability to combine those samples:

<iframe
  src="./assets/bernstein-convergence.html"
  title="Bernstein polynomial convergence at x equals 0.75"
  width="100%"
  height="820"
  loading="lazy"
  sandbox="allow-scripts"
></iframe>

[Open the Bernstein convergence explorer](./assets/bernstein-convergence.html) · [Open the published version](https://bernstein-convergence.priyanshu-arora2007.chatgpt.site/)

## The theorem

**Weierstrass approximation theorem.** Let $f:[0,1]\to\mathbb{R}$ be continuous. For every $\varepsilon>0$, there is a polynomial $p$ such that

\[
\sup_{x\in[0,1]} |p(x)-f(x)| < \varepsilon.
\]

Equivalently, there is a sequence of polynomials $p_n$ that converges uniformly to $f$ on $[0,1]$.

## The Bernstein polynomials

Divide $[0,1]$ into the grid

\[
0,\frac1n,\frac2n,\ldots,\frac{n-1}{n},1
\]

and sample $f$ at those points. Define

\[
B_n f(x)
=
\sum_{k=0}^{n}
f\!\left(\frac{k}{n}\right)
\binom{n}{k}x^k(1-x)^{n-k}.
\]

For each fixed $n$, this is a polynomial in $x$. We will prove that

\[
B_n f \longrightarrow f
\qquad\text{uniformly on }[0,1].
\]

The grid values $f(k/n)$ act as anchors, but $B_n f$ is not generally an interpolation polynomial: it does not have to pass through every interior anchor. Instead, it forms a weighted average of all the anchor heights.

## The probabilistic interpretation

Fix $x\in[0,1]$. Flip a coin $n$ times, with probability $x$ of heads on each flip. Let $X_1,\ldots,X_n$ be the indicator variables for the individual flips, and let

\[
S_n=X_1+\cdots+X_n.
\]

Then $S_n\sim\operatorname{Binomial}(n,x)$, so

\[
\Pr(S_n=k)=\binom nk x^k(1-x)^{n-k}.
\]

Consequently,

\[
B_n f(x)
=
\mathbb E\!\left[f\!\left(\frac{S_n}{n}\right)\right].
\]

The random point $S_n/n$ has mean $x$:

\[
\mathbb E\!\left[\frac{S_n}{n}\right]=x.
\]

Its variance is

\[
\operatorname{Var}\!\left(\frac{S_n}{n}\right)
=
\frac{x(1-x)}{n}
\leq
\frac{1}{4n}.
\]

Thus $S_n/n$ becomes increasingly concentrated near $x$ as $n$ grows. The value $B_n f(x)$ is therefore an average of values of $f$ taken mostly at grid points close to $x$.

## Proof of uniform convergence

Let $\varepsilon>0$.

Because $f$ is continuous on the compact interval $[0,1]$, it is uniformly continuous. Therefore there is a number $\delta>0$ such that, for every $s,t\in[0,1]$,

\[
|s-t|<\delta
\quad\Longrightarrow\quad
|f(s)-f(t)|<\frac{\varepsilon}{2}.
\]

The same $\delta$ works everywhere on the interval.

Also, because $f$ is continuous on a compact interval, it is bounded. Choose $M\geq 0$ such that

\[
|f(t)|\leq M
\qquad\text{for every }t\in[0,1].
\]

For an arbitrary $x\in[0,1]$, use the probabilistic representation to write

\[
\begin{aligned}
|B_n f(x)-f(x)|
&=
\left|
\mathbb E\!\left[
f\!\left(\frac{S_n}{n}\right)-f(x)
\right]
\right| \\
&\leq
\mathbb E\!\left[
\left|f\!\left(\frac{S_n}{n}\right)-f(x)\right|
\right].
\end{aligned}
\]

Now a major theme of probability is partitions to divide the sample space into what we (Mathematicians) want. Split the expectation according to whether the random point $S_n/n$ is close to $x$. Define

\[
G=\left\{\left|\frac{S_n}{n}-x\right|<\delta\right\},
\qquad
G^c=\left\{\left|\frac{S_n}{n}-x\right|\geq\delta\right\}.
\]

Then

\[
\begin{aligned}
|B_n f(x)-f(x)|
&\leq
\mathbb E\!\left[
\left|f\!\left(\frac{S_n}{n}\right)-f(x)\right|
\mathbf 1_G
\right] \\
&\quad+
\mathbb E\!\left[
\left|f\!\left(\frac{S_n}{n}\right)-f(x)\right|
\mathbf 1_{G^c}
\right].
\end{aligned}
\]

### The good event

On $G$, the point $S_n/n$ lies within $\delta$ of $x$. Uniform continuity gives

\[
\left|f\!\left(\frac{S_n}{n}\right)-f(x)\right|
<\frac{\varepsilon}{2}.
\]

Therefore

\[
\mathbb E\!\left[
\left|f\!\left(\frac{S_n}{n}\right)-f(x)\right|
\mathbf 1_G
\right]
\leq \frac{\varepsilon}{2}.
\]

### The bad event

On $G^c$, the two function values might be far apart, but boundedness gives the crude estimate

\[
\left|f\!\left(\frac{S_n}{n}\right)-f(x)\right|
\leq 2M.
\]

Hence

\[
\mathbb E\!\left[
\left|f\!\left(\frac{S_n}{n}\right)-f(x)\right|
\mathbf 1_{G^c}
\right]
\leq
2M\Pr(G^c).
\]

Chebyshev's inequality gives

\[
\begin{aligned}
\Pr(G^c)
&=
\Pr\!\left(
\left|\frac{S_n}{n}-x\right|\geq\delta
\right) \\
&\leq
\frac{\operatorname{Var}(S_n/n)}{\delta^2} \\
&=
\frac{x(1-x)}{n\delta^2} \\
&\leq
\frac{1}{4n\delta^2}.
\end{aligned}
\]

It follows that the contribution from the bad event is at most

\[
2M\Pr(G^c)
\leq
\frac{M}{2n\delta^2}.
\]

### Combining the estimates

For every $x\in[0,1]$,

\[
|B_n f(x)-f(x)|
\leq
\frac{\varepsilon}{2}
+
\frac{M}{2n\delta^2}.
\]

Choose an integer $N$ large enough that

\[
N>\frac{M}{\varepsilon\delta^2}.
\]

Then, whenever $n\geq N$,

\[
\frac{M}{2n\delta^2}<\frac{\varepsilon}{2},
\]

and therefore

\[
|B_n f(x)-f(x)|<\varepsilon
\qquad\text{for every }x\in[0,1].
\]

Taking the supremum over $x$ gives

\[
\sup_{x\in[0,1]}|B_n f(x)-f(x)|<\varepsilon.
\]

This proves that $B_n f\to f$ uniformly on $[0,1]$. □

## Why fixing $x$ still proves uniform convergence

The probabilistic interpretation begins by fixing $x$, but the final choice of $N$ does not depend on $x$. That independence comes from three uniform facts:

1. Uniform continuity provides one $\delta$ that works for every point of $[0,1]$.
2. Boundedness provides one $M$ that bounds $f$ everywhere.
3. The variance estimate uses $x(1-x)\leq 1/4$, which holds for every $x\in[0,1]$.

The order of the quantifiers is therefore

\[
\forall\varepsilon>0\;\exists N\;\forall n\geq N\;\forall x\in[0,1],
\qquad
|B_n f(x)-f(x)|<\varepsilon.
\]

If the required $N$ depended on $x$, the argument would prove only pointwise convergence. Here it does not, which is precisely why the convergence is uniform.

## Visual intuition

We basically chopped the interval [0,1 ] into n slices

At a chosen point $x$, the coefficients

\[
\binom nk x^k(1-x)^{n-k}
\]

form a probability distribution over the grid points $k/n$. The Bernstein polynomial averages the anchor heights $f(k/n)$ using these weights.

As $n$ increases:

- the grid becomes finer;
- the weights concentrate more tightly around $k/n\approx x$;
- nearby anchors have heights close to $f(x)$, by continuity; and
- faraway anchors may have very different heights, but their total weight tends to zero.

The Bernstein polynomial is therefore a smoothly blurred version of the sampled function, with the blur shrinking uniformly as $n\to\infty$.
