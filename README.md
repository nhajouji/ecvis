# Elliptic Curves over Finite Fields: Visualizations and Data Science

## Introduction

We're going to play around with two seemingly unrelated types of mathematical objects:

1. Elliptic curves over finite fields

2. Lattices in the complex plane

## Key Math Ideas

Let $E/\mathbb{F}_p$ be an elliptic curve and let $\phi_E : E \to E$ be the Frobenius $\phi: (x,y)\mapsto (x^p,y^p)$.
Our goal is to visualize this curve.

Visualizing things in positive characteristic is always a bit weird, because  pictures live in $\mathbb{R}^n$ and you simply can't fit a finite field in a real vector space in any meaningful way.
However, with elliptic curves, enough crazy things work out that we can actually produce meaningful pictures by lifting everything (in a meaningful way) to characteristic 0.
The magic result that makes this possible is the following: we can lift the pair $(E, \phi)$ to a pair $(\Lambda, \alpha)$,
where:
* $\mathbb{C}/\Lambda$ is isomorphic to an algebraic elliptic curve $E^{alg}/\mathbb{C}$,
where $E^{alg}$ is given by an equation with algebraic integer coefficients.
* The reduction of $E^{alg}$ mod $p$ is isomorphic to $E$ as elliptic curves over $\mathbb{F}_p$.
* The endomorphism $z + \Lambda \mapsto \alpha z + \Lambda$ reduces to $\phi_E$ mod $p$.

Once we know the lift $(\Lambda, \alpha)$,
we're construct the superlattice $\Lambda_\alpha = \frac{1}{\alpha-1}\Lambda \supseteq \Lambda$.
The quotient $\Lambda_\alpha/\Lambda \subset \mathbb{C}/\Lambda$, in the algebraic world, corresponds a finite subgroup of $E^{alg}$.

* If $E^{alg}$ is an elliptic curve over the algebraic integers,
and $P = (x,y)$ is a point of finite order on $E$,
then $x,y$ must also be algebraic integers.

* Let $E^{alg}(\overline{\mathbb{Z}})$ be the set of points on $E^{alg}$ with algebraic integer coordinates.
The previous bullet says that $E^{alg}(\overline{\mathbb{Z}})$ contains all points of finite order, 
and in particular, contains the subset $\Lambda_\alpha/\Lambda$.

* Let $p$ be a prime and let $\mathfrak{m} \subset \overline{\mathbb{Z}}$ be a maximal ideal containing $p$.
The quotient $\overline{\mathbb{Z}}/\mathfrak{m}$ is isomorphic to the algebraic closure of $\mathbb{F}_p$.
The reduction mod $p$ map is the map:

  $ E^{alg}(\overline{\mathbb{Z}}) \to E^{alg}(\overline{\mathbb{Z}}/\mathfrak{m}) \cong E ( \overline{\mathbb{F}_p})$
  
that takes $(x,y)$ to $(x\pmod \mathfrak{m},y\pmod\mathfrak{m})$.

We will mainly be interested in the restriction of this map to $E^{alg}_{tors}$.


* If the model $E^{alg}$ reduces to a smooth elliptic curve mod $p$,
then the reduction map, restricted to the subgroup of points of finite order
$E^{alg}_{tors} \to E(\overline{\mathbb{F}_p})$,
is injective.
In particular, the map $\Lambda_\alpha/\Lambda \to E \pmod p$ is injective.

* Let $n$ be an integer and $E(\mathbb{F}_{p^n})$ the subgroup of points defined over $\mathbb{F}_{p^n}$.
Then $E(\mathbb{F}_{p^n})$ coincides with the image of $\Lambda_{\alpha^n}/\Lambda$ under the reduction map.
(And since the reduction map is injective, those groups are isomorphic)

The upshot is this: we can easily make a picture that includes as many $\Lambda_n$ as we wish,
and we can interpret that picture as depicting the points on $E$ defined over the corresponding extensions of $\mathbb{F}_p$.

