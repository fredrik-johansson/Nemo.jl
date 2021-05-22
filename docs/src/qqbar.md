```@meta
CurrentModule = Nemo
```

# Algebraic numbers

Nemo allows working with exact real and complex algebraic numbers.

The default algebraic number type in Nemo is provided by Calcium. The
associated field of algebraic numbers is represented by the constant
parent object called `CalcumQQBar`.

For convenience we define

```
QQBar = CalciumQQBar
```

so that algebraic numbers can be constructed using `QQBar` instead of
`CalciumQQBar`. Note that this is the name of a specific parent object,
not the name of its type.


 Library        | Element type  | Parent type
----------------|---------------|--------------------
Calcium         | `qqbar`       | `CalciumQQBarField`

## Important note on performance

The default algebraic number type represents algebraic numbers
in canonical form using minimal polynomials. This works well for representing
individual algebraic numbers, but it does not provide the best
performance for field arithmetic.
To compute in a fixed subfield of $\overline{\mathbb{Q}}$,
it is typically far more efficient to fix a generator $a$
and construct an Antic number field to represent $\mathbb{Q}(a)$.

In the future, Nemo will provide an alternative implementation
of algebraic numbers using Calcium field elements instead of minimal
polynomials. The minimal polynomial representation will still be available
since it provides canonical forms with predictable behavior.

## Algebraic number functionality

### Examples

Solving the quintic equation:

```julia
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rational Field, x)

julia> v = roots(x^5-x-1, CalciumQQBar)
5-element Array{qqbar,1}:
 Root 1.16730 of x^5 - x - 1
 Root 0.181232 + 1.08395*I of x^5 - x - 1
 Root 0.181232 - 1.08395*I of x^5 - x - 1
 Root -0.764884 + 0.352472*I of x^5 - x - 1
 Root -0.764884 - 0.352472*I of x^5 - x - 1

julia> v[1]^5 - v[1] - 1 == 0
true
```

