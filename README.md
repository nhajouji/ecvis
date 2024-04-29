# Elliptic Curves over Finite Fields: Visualizations and Data Science

## Introduction

We're going to play around with two seemingly unrelated types of mathematical objects:

1. Elliptic curves over finite fields

2. Lattices in the complex plane

## Key Math Ideas

Let $E/\mathbb{F}_p$ be an elliptic curve and let $\phi : E \to E$ be the Frobenius $\phi: (x,y)\mapsto (x^p,y^p)$.
Our goal is to visualize this curve.

Visualizing things in positive characteristic is always a bit weird, because  pictures live in $\mathbb{R}^n$ and you simply can't fit a finite field in a real vector space in any meaningful way.
However, with elliptic curves, enough crazy things work out that we can actually produce meaningful pictures by lifting everything (in a meaningful way) to characteristic 0.

The magic result that makes this possible is the following: We can always find a curve $\tilde{E}$, together with an endomrphism $\tilde{\phi}: \tilde{E}\to \tilde{E}$, such that the reduction of $\tilde{E}$ mod $p$ is $E$ and the reduction of $\tilde{\phi}$ is the Frobenius.

We will use this to lift the points of $E$ that are defined over $\mathbb{F}_p$ (or more generally the points that are defined over $\mathbb{F}_{p^n}$) to points on $\tilde{E}$.

* The curve $\tilde{E}$ is necessarily isogenous to a curve with CM so its $j$-invariant is an algebraic integer. Thus, we can find a model of $\tilde{E}$ which is defined over the algebraic integers, and thus can be reduced modulo any prime. 
* Nextt, note that every point which is defined over a finite field necessarily has finite order in the Mordell-Weil group - thus, we only need to look for torsion points on $\tilde{E}$ that reduce to points which are defined over our finite field of choice ($\mathbb{F}_p$ or $\mathbb{F}_{p^m}$).
*  The reduction map $\tilde{E}(\overline{\mathbb{Q}})_{tors} \to \tilde{E}(\overline{\mathbb{F}_p})$ is surjective and almost injective - the only points in the kernel are points of order $p^m$, because elliptic curves in characteristic $p$ either have no points of order $p$ or the $p^\infty$ subgroup has rank 1. 
* The key thing is that every point defined over a finite field gives rise to an essentially unique torsion point on $\tilde{E}$ which has the same order in the Mordell-Weil group. The only possible ambiguity that could arise would be if there is a point of order $p$ in characteristic $p$ (although this won't be an issue).
* To determine whether a given torsion point $\tilde{P}$ on $\tilde{E}$ reduces, in characteristic $p$, to a point defined over $\mathbb{F}_{p^m}$, all we have to do is check whether $\tilde{\phi}^m(\tilde{P}) = \tilde{P}$.
* In fact, we don't need to go point by point - we can simply compute the kernel of $\tilde{\phi}^m - \mathrm{id}$. These are exactly the points on $\tilde{E}$ that are fixed by $\tilde{\phi}^m$ and correspond to the points in characteristic $p$ that are fixed by $(x,y)\mapsto (x^{p^m},y^{p^m})$.

To obtain our visualizations, we will use lifts of $E$ of the form $\mathbb{C}/\Lambda$. On an analytic model of this form, the lift of Frobenius is just a map of the form $z+\Lambda \mapsto \alpha z + \Lambda$, where $\alpha$ is a complex algebraic integer of norm $p$.

The set of points on $\CC/\Lambda$ that are fixed by $\tilde{\phi}^m : z+\Lambda \mapsto \alpha^m z + \Lambda$ is simply $\frac{1}{\alpha^m-1}\Lambda + \Lambda$, so once we know $\alpha, m,\Lambda$ we can immediately obtain pictures.





## Lattices in the complex plane

A lattice in the plane is just a collection of regularly spaced points that go on forever.
Geometrically, they are like higher dimensional analogs of the set of integers on the number line.

Algebraically, we can think of lattices as "discrete vector spaces":
* Every lattice can be described using a basis, which will be some linearly independent subset of vectors in a continuous vector space.
* The associated lattice consists of the set of points that can be obtained as a linear combination of the basis elements, using only INTEGER linear combinations.
  
A lattice in the complex plane is just a subset of points that can be obtained as integer linear combinations of some pair of complex numbers $\tau_1, \tau_2$.

### Doubly periodic functions

Let $\Lambda$ be a lattice in $\mathbb{C}$.
A $\Lambda$-periodic function is a function $f : \mathbb{C} \to \mathbb{C} \cup \{ \infty\}$ 
with the property that $f(z+\lambda) = f(z)$ for all $z \in \mathbb{C}$ and all $\lambda \in \Lambda$.

The most important example of a $\Lambda$-periodic function is the 'Weierstrass $\wp$-function':

$$ \wp_\Lambda(z) = \frac{1}{z^2} + \sum_{\lambda \in \Lambda \\ lambda \neq 0} \frac{1}{(z-\lambda)^2}-\frac{1}{\lambda^2}$$

The Weierstrass $\wp$-function plays the role that sine/cosine play in the theory of "normal" periodic functions.

Just as sine, cosine satisfy a fundamental relation - $\sin^2 +\cos^2 = 1$ - the Weierstrass function and its derivative satisfy an equation of the form:
$$\wp_\Lambda'(z)^2 = 4\wp_\Lambda(z)^3 + f(\Lambda) \wp_\Lambda(z) + g(\Lambda) $$
where $f(\Lambda), g(\Lambda)$ are fixed complex numbers that depend on $\Lambda$.


### The $j$-map
Two lattices are "homothetic" if one can be transformed into the other using rotations and rescalings. For a pair of lattices $\Lambda_1, \Lambda_2$ in the complex plane, this is equivalent to saying $\Lambda_1 = z \Lambda_2$ for some complex number $z$.
(Here $z \Lambda_2$ means the set obtained by multiplying every element in the lattice $\Lambda_2$ by the complex number $z$).

Given $\tau_1, \tau_2$, we can obtain a new basis that represents a homothetic lattice by dividing by one of the two complex numbers - i.e. we can take $\frac{\tau_1}{\tau_1} = 1, \frac{\tau_2}{\tau_1}$ or $\frac{\tau_2}{\tau_1}, \frac{\tau_2}{\tau_2} = 1$. This is convenient, since one of the new basis vectors will be 1, so we only have to keep track of the second. Thus, it is common to assume that the basis is $1,\tau$. Furthermore, swapping the roles of $\tau_1, \tau_2$ if necessary, we may assume $\tau$ is an element of the upper half plane.

We will view points on the upper half plane as representing lattices in the complex plane.

Even though we've eliminated a few degrees of freedom, there are still many elements $\tau$ we can use to represent any lattice. For example, replacing $\tau$ by $\tau' = \tau+1$ does not change the lattice, because any integral linear combination we can obtain using $1, \tau$ can be obtained using $1,\tau+1$ and vice versa.

Two elements $\tau_1, \tau_2$ represent the same homothety class of lattices if and only if we can find $a,b, c, d$ satisfying:
* $ad-bc = 1$
* $\frac{a\tau_1+b}{c\tau_1+d} = \tau_2$.

To decide whether $\tau_1, \tau_2$ represent homothetic lattices, we have the following tools:
* There are practical algorithms that allow us to find such an $a,b,c,d$, if any exist.
* There is also a surjective function $j_{an}: \mathcal{H}\to \mathbb{C}$ with the following property:
$j_{an}(\tau) = j_{an}(\tau')$ if and only if the lattices generated by $1,\tau$ and $1,\tau'$ are homothetic.

Now, we can't actually compute the function $j_{an}$ explicitly, but we will be able to do something almost as good:
For (infinitely many, but only countably infinitely many) lattices, the value of $j_{an}(\tau)$ will be an "algebraic integer". This means $j_{an}(\tau)$ can be described as a root of a monic polynomial with integer coefficients.
Once we obtain the polynomial associated to $\tau$, we can solve it to obtain the value of $j_{an}(\tau)$.

There's only one problem: the polynomial will have multiple roots. Each of these other roots comes from a lattice generated by some $\tau_i$ which is related to the original $\tau$.

## Elliptic Curves

An elliptic curve is the set of solutions to an equation of the form:
$$ y^2 = x^3 + fx + g$$

For example, the image of $\mathbb{C}$ under the map $z \mapsto (\wp(z),\wp'(z))$ is an elliptic curve.




## Elliptic Curves over Finite Fields

