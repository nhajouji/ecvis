# Visualizing Elliptic Curves (and other algebraic groups) over Finite Fields

## Introduction

Our goal is to make meaningful illustrations of varieties over finite fields.

For concreteness, suppose we have a plane curve in $\mathbb{F}_p^2$. There is an obvious way of making a picture that represents the curve:
* We represent the ambient space as a $p \times p$ grid of points.
* We color in the points on the variety in a different color.

These pictures don't shed any insight on the curve - we may as well just list the points.

A different way to approach this problem is by trying to "lift" the variety over $\mathbb{F}_p$ to a subset of a variety in characteristic 0, assuming we have some way of visualizing the characteristic 0 variety.
This can be done, but without adding additional constraints, the picture will not be any more meaningful than the previous one - there are infinitely many ways of lifting something in characteristic $p$ to characteristic 0 and they all look different.
However, this idea can lead to meaningful illustrations if try to lift not just the points of the variety, but also some additional structure.

## Commutative algebraic groups over finite fields
### Reduction
Let:
* $p \in \mathbb{Z}$ is a prime
* $\mathbb{Z}^{alg}$ is the ring of algebraic integers.

In order to define lifts, we first need to extend the definition of ``mod $p$" to the ring of algebraic integers.

This is straightforward: we just need to take a quotient by a maximal ideal of $\mathbb{Z}^{alg}$ that contains the prime $p$.
The quotient ring is isomorphic to the algebraic closure of $\mathbb{F}_p$,
so we can think of the quotient map as a surjection of rings
 $\rho : \mathbb{Z}^{alg} \to \mathbb{F}_p^{alg}$.
 We say that $\xi_1 \equiv \xi_2 \pmod p$ if $\rho(\xi_1) = \rho(\xi_2)$.

However, it is important to note that the map $\rho$ is not unique - we will be forced to deal with that non-uniqueness when we work with elliptic curves.
It might be better to write $\pmod \rho$.


### Lifting
Let $X/\mathbb{F}_p$ be a commutative algebraic group,
and let $G$ be the Galois group $Gal(\mathbb{F}_p^{sep}/\mathbb{F}_p)$.
* The group $G$ is generated by the Frobenius automorphism $x\mapsto x^p$.
* $G$ acts on $X$ by endomorphisms of $X$ a a commutative algebraic group (i.e. the maps $x\mapsto g \cdot x$ are morphisms of $X$ as a group and a variety.)
* To describe the action of $G$, we just need to know the endomorphism $\phi : X \to X$ induced by the Frobenius. Note that $\phi$ fixes the points of $X$ that are defined over $\mathbb{F}_p$, and permutes the points of $X$ that are defined over $\mathbb{F}_p^{sep}$.

Given $X$ as above, and a reduction map $\rho$,
we define an algebraic lift of $X$ to be a pair $(\tilde{X}, \tilde{\phi})$,
where:
* $\tilde{X}$ is a commutative algebraic group whose reduction mod $\rho$ is isomorphic to $X$.
* $\tilde{\phi}: \tilde{X}\to \tilde{X}$ is a morphism of commutative algebraic groups that reduces to $\phi$ mod $\rho$.

### What we can hope to visualize using lifts
Given a lift $\tilde{X}, \tilde{\phi}$,
define $\tilde{X}(\mathbb{F}_{p^n})$ to be the set of fixed points of
$\tilde{\phi}^n$.
Then:
* The group homorphism $\tilde{X}(\mathbb{Z}^{alg})\to X(\mathbb{F}_p^{sep})$ restricts to an isomorphism $\tilde{X}(\mathbb{F}_{p^n})\to X(\mathbb{F}_{p^n})$,
so we can use pictures of $\tilde{X}(\mathbb{F}_{p^n})$ to visualize $X(\mathbb{F}_{p^n})$. 
* $\tilde{X}(\mathbb{F}_{p^n})$ is a subgroup (it coincides with the kernel of $\tilde{\phi}^n - \mathrm{id}$.), so this allows us to visualize the group structure better (assuming we can visualize it in characteristic 0).
* As a bonus, we can also see how Frobenius acts on the points of $\tilde{X}$ that are defined over any finite extension of $\mathbb{F}_p$.


#### Uniqueness

In a nutshell, here is why the pictures we obtain are "unique":
* The groups $X(\mathbb{F}_{p^n})$ are all finite, so every point has finite order.
* In characteristic 0, we will be able to use arguments like 'This is the only subgroup of order 6 so it must be the correct subgroup'.

I will be more precise when we get to examples (next).

## Examples

We're going to focus on 1-dimensional algebraic groups: $X$ is either the additive group $\mathbb{G}_a$,
the multiplicative group $\mathbb{G}_m$ or an elliptic curve.
* We're not actually going to try picturing the additive group $\mathbb{G}_a$, because it does not have any algebraic lifts: this is because every element in the additive group in characteristic $p$ has finite order, but the additive group is torsionfree in characteristic 0. I'm only mentioning it to highlight that when lifts *do* exist, something special is happening!

* Lifts of the multiplicative group *do* exist, and they are super easy to understand! The multiplicative group of a finite field is cyclic, and the multiplicative group in characteristic 0 has exactly one cyclic subgroup of order $n$ for each $n \in \mathbb{N}$ (generated by a primitive $n$th root of unity),
so the whole story is completely straightforward. The lift of Frobenius to characteristic 0 is given by the same formula ($x\mapsto x^p$).

Lifts of elliptic curves are where things get interesting:
* They always exist.
* We can use make pictures using very specific lattices that vary from curve to curve.
* Assuming a reduction map has been fixed, lifts are unique up to isomorphism.
* *But* there is a weird ambiguity - changing the reduction map can permute which lift belongs to which elliptic curve $\pmod p$.

We will explain this "ambiguity" next.

## Algebraic lattices and Analytic Lifts

Let $X(1) = \mathcal{H}/SL_2(\mathbb{Z})$ be the quotient of the upper half plane by the modular group, and $j :  X(1) \to \mathbb{CP}^1$ the classical $j$-invariant.

Let $\Lambda \subset \mathbb{C}$ be a lattice - the homothety class of $\Lambda$ gives rise to a point on $X(1)$,
and the $j$-invariant sends the homothety class of $\Lambda$ to an element of $\mathbb{CP}^1$.

* We say that $\Lambda$ is algebraic if $j([\Lambda])$ is an algebraic number.
* We say that $\Lambda$ is integral if $j([\Lambda])$ is an algebraic integer.
* If $\Lambda, \Lambda'$ are algebraic lattices, we say they are algebraically similar if their $j$-invariants have the same minimal polynomial over $\mathbb{Z}$,
i.e. if the $j$-invariants are Galois conjugates of each other.

### Analytic lifts

Fix a prime $p$ and a reduction map $\rho : \mathbb{Z}^{alg} \to \mathbb{F}_p^{alg}$, and let $E/\mathbb{F}_p$ be an elliptic curve.

An analytic lift of $E$ is a pair $(\Lambda, \alpha)$,
where:
* $\Lambda$ is an algebraic integral lattice, 
* $\alpha$ is a complex number satisfying $\alpha \Lambda \subset \Lambda$, and 
* The pair $(\mathbb{C}/\Lambda, z+\Lambda\mapsto \alpha z + \Lambda)$ (which consists of an elliptic curve over the algebraic integers, and an endomorphism of that elliptic curve) is an algebraic lift of $E$ in the sense we defined earlier. 

We say that two elliptic curves $E, E'/\mathbb{F}_p$ are analytically similar if the following is true: there is a pair $(\Lambda,\alpha)$, and a pair of reduction maps $\rho, \rho'$, such that $(\Lambda, \alpha)$ is a lift of $E$ wrt $\rho$,
and a lift of $E'$ wrt $\rho'$.

### Analytic vs Algebraic similarity

Fix a reduction map $\rho_0$, and let $E, E'$ be elliptic curves.

Suppose $E, E'$ are analytically similar, and let $(\Lambda, \alpha), (\Lambda, \alpha')$ be analytic lifts wrt to $\rho_0$. Then the lattices $\Lambda, \Lambda'$ must be algebraically similar (and $\alpha = \alpha'$).

On the other hand, suppose $(\Lambda, \alpha)$ is a lift of $E$ wrt to $\rho_0$, and let $\Lambda'$ be a lattice which is algebraically equivalent to $\Lambda$. Then:
* The complex number $\alpha$ satisfies $\alpha \Lambda' \subset \Lambda'$.
* The pair $(\Lambda', \alpha)$ is an analytic lift of some $E'/\mathbb{F}_p$ which is analytically similar to $E$.

### Takeaway

It does not make sense to say "This picture represents this specific curve" without specifying the reduction map.

However, we can say "This algebraic equivalence class of lattices represents this analytic equivalence class of elliptic curves" without needing to specify the reduction map. (Also, it is possible that the equivalence class contains a single element, and if that's the case, we can say "this picture represents this curve")

So, when we make pictures, we will be making several pictures at a time - and we can't really say which of the pictures represents which of the polynomial equations. 

