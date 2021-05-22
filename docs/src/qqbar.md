```@meta
CurrentModule = Nemo
```

# Algebraic numbers

Nemo allows working with exact real and complex algebraic numbers.

The default algebraic number type in Nemo is provided by Calcium. The
associated field of algebraic numbers is represented by the constant
parent object called `CalciumQQBar`.

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

**Important note on performance**

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

### Constructing algebraic numbers

```@docs
rand(R::CalciumQQBarField; degree::Int, bits::Int, randtype::Symbol=:null)
```

### Conversions

**Examples**

Integer and rational algebraic numbers can be converted to Nemo
integer and rational types. Algebraic numbers can be evaluated
numerically to arbitrary precision by converting
to real or complex Arb fields:

```julia
julia> fmpz(QQBar(3))
3

julia> fmpq(QQBar(3) // 2)
3//2

julia> RR = ArbField(64); RR(sqrt(QQBar(2)))
[1.414213562373095049 +/- 3.45e-19]

julia> CC = AcbField(32); CC(QQBar(-1) ^ (QQBar(1) // 4))
[0.707106781 +/- 2.74e-10] + [0.707106781 +/- 2.74e-10]*im
```

### Minimal polynomials and conjugates

**Examples**

Solving the quintic equation:

```julia
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rational Field, x)

julia> v = roots(x^5-x-1, QQBar)
5-element Array{qqbar,1}:
 Root 1.16730 of x^5 - x - 1
 Root 0.181232 + 1.08395*I of x^5 - x - 1
 Root 0.181232 - 1.08395*I of x^5 - x - 1
 Root -0.764884 + 0.352472*I of x^5 - x - 1
 Root -0.764884 - 0.352472*I of x^5 - x - 1

julia> v[1]^5 - v[1] - 1 == 0
true
```

Computing exact eigenvalues of a matrix:

```julia
julia> eigenvalues(ZZ[1 1 0; 0 1 1; 1 0 1], QQBar)
3-element Array{qqbar,1}:
 Root 2.00000 of x - 2
 Root 0.500000 + 0.866025*I of x^2 - x + 1
 Root 0.500000 - 0.866025*I of x^2 - x + 1
```

Retrieving the minimal polynomial and algebraic conjugates
of a given algebraic number:

```julia
julia> minpoly(PolynomialRing(ZZ, "x")[1], QQBar(1+2im))
x^2 - 2*x + 5

julia> conjugates(QQBar(1+2im))
2-element Array{qqbar,1}:
 Root 1.00000 + 2.00000*I of x^2 - 2x + 5
 Root 1.00000 - 2.00000*I of x^2 - 2x + 5
```

**Interface**

```@docs
minpoly(R::FmpzPolyRing, x::qqbar)
minpoly(R::FmpqPolyRing, x::qqbar)
roots(f::fmpz_poly, R::CalciumQQBarField)
roots(f::fmpq_poly, R::CalciumQQBarField)
conjugates(a::qqbar)
eigenvalues(A::fmpz_mat, R::CalciumQQBarField)
eigenvalues(A::fmpq_mat, R::CalciumQQBarField)
```

### Comparing algebraic numbers

### Roots of unity and trigonometric functions

**Examples**

```julia
julia> sinpi(QQBar(7) // 13)
Root 0.992709 of 4096x^12 - 13312x^10 + 16640x^8 - 9984x^6 + 2912x^4 - 364x^2 + 13

julia> tanpi(atanpi(sqrt(QQBar(2)) + 1))
Root 2.41421 of x^2 - 2x - 1

julia> root_of_unity(QQBar, 5)
Root 0.309017 + 0.951057*I of x^4 + x^3 + x^2 + x + 1

julia> root_of_unity(QQBar, 5, 4)
Root 0.309017 - 0.951057*I of x^4 + x^3 + x^2 + x + 1

julia> w = (1 - sqrt(QQBar(-3)))//2
Root 0.500000 - 0.866025*I of x^2 - x + 1

julia> is_root_of_unity(w)
true

julia> is_root_of_unity(w + 1)
false

julia> root_of_unity_as_args(w)
(6, 5)
```

**Interface**

```@docs
root_of_unity(C::CalciumQQBarField, n::Int)
root_of_unity(C::CalciumQQBarField, n::Int, k::Int)
is_root_of_unity(a::qqbar)
root_of_unity_as_args(a::qqbar)
exp_pi_i(a::qqbar)
log_pi_i(a::qqbar)
sinpi(a::qqbar)
cospi(a::qqbar)
tanpi(a::qqbar)
asinpi(a::qqbar)
acospi(a::qqbar)
atanpi(a::qqbar)
```

### Guessing

**Examples**

An algebraic number can be recovered from a numerical value:

```julia
julia> RR = ArbField(53); guess(QQBar, RR("1.41421356 +/- 1e-6"), 2)
Root 1.41421 of x^2 - 2
```

**Interface**

```@docs
guess(R::CalciumQQBarField, x::arb, maxdeg::Int, maxbits::Int=0)
guess(R::CalciumQQBarField, x::acb, maxdeg::Int, maxbits::Int=0)
```

