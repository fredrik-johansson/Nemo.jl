```@meta
CurrentModule = Nemo
```

# Exact real and complex numbers

Exact real and complex numbers are provided by Calcium.
Internally, a number $z$ is represented as an element
of an extension field of the rational numbers. That is,

$z \in \mathbb{Q}(a_1,\ldots,a_n)$

where $a_1, \ldots, a_n$
are symbolically defined algebraic or transcendental
real or complex numbers
such as $\pi$, $\sqrt{2}$ or $e^{\sqrt{2} \pi i}$.
The user does not normally need to worry about the details of
the internal representation;
Calcium constructs extension numbers and fields automatically
as needed to perform operations.

The user must create a `CalciumField` instance which represents the
mathematical domain $\mathbb{C}$. This parent object holds a cache of extension
numbers and fields used to represent individual elements. It also
stores various options for evaluation (documented further below).


 Library        | Element type  | Parent type
----------------|---------------|--------------------
Calcium         | `ca`          | `CalciumField`


Please note the following:

* It is in the nature of exact complex arithmetic that some operations
  must be implemented using incomplete heuristics. For example, testing
  whether an element is zero will not always succeed. When Calcium is
  unable to perform a task, Nemo will throw an exception.
  This ensures that Calcium fields behave exactly
  and never silently return wrong results.

* Calcium elements can optionally hold special non-numerical values:

  * Unsigned infinity $\hat \infty$

  * Signed infinities ($\pm \infty$, $\pm i \infty$, and more generally $e^{i \theta} \cdot \infty$)

  * Undefined

  * Unknown

  By default, such special values are
  disallowed so that a `CalciumField` represents the
  mathematical field $\mathbb{C}$, and any operation that would result
  in a special value (for example, $1 / 0 = \hat \infty$) will throw an exception.
  To allow special values, pass `extended=true` to the `CalciumField` constructor.

* `CalciumField` instances only support single-threaded use.
  You must create a separate parent object for each thread
  to do parallel computation.

* When performing an operation involving two `ca` operands with different
  parent objects, Nemo will arbitrarily coerce the operands (and hence
  the result) to one of the parents.



## Calcium field options

The `CalciumField` parent stores various options that affect
simplification power, performance, or appearance.
The user can override any of the default values using
`C = CalciumField(options=dict)` where `dict` is a dictionary
with `Symbol => Int` pairs. To retrieve the option values as a dictionary
(including any default values not set by the user),
call `options(C)`.

The following options are supported:

 Option                   | Explanation
--------------------------|------------------------------------------------
  `:verbose`              | Enable debug output
  `:print_flags`          | Flags controlling print style
  `:mpoly_ord`            | Monomial order for polynomials
  `:prec_limit`           | Precision limit for numerical evaluation
  `:qqbar_deg_limit`      | Degree limit for algebraic numbers
  `:low_prec`             | Initial precision for numerical evaluation
  `:smooth_limit`         | Factor size limit for smooth integer factorization
  `:lll_prec`             | Precision for integer relation detection
  `:pow_limit`            | Maximum exponent for in-field powering
  `:use_gb`               | Enable Gröbner basis computation
  `:gb_length_limit`      | Maximum ideal basis length during Gröbner basis computation
  `:gb_poly_length_limit` | Maximum polynomial length during Gröbner basis computation
  `:gb_poly_bits_limit`   | Maximum bit size during Gröbner basis computation
  `:gb_vieta_limit`       | Maximum degree to use Vieta's formulas
  `:trig_form`            | Default form of trigonometric functions

An important function of these options is to control how hard Calcium will
try to find an answer before it gives up. For example:

* Setting `:prec_limit => 65536` will allow Calcium to use up to 65536
  bits of precision (instead of the default 4096) to prove inequalities.

* Setting `:qqbar_deg_limit => typemax(Int)` (instead of the default 120)
  will force most calculations involving algebraic numbers to run to
  completion, no matter how long this will take.

* Setting `:use_gb => 0` (instead of the default 1) disables use of
  Gröbner bases. In general, this will negatively impact Calcium's ability
  to simplify field elements and prove equalities, but it can speed up
  calculations where Gröbner bases are unnecessary.

For a detailed explanation, refer to the following section
in the Calcium documentation:
<https://fredrikj.net/calcium/ca.html#context-options>

## Basic examples

```julia
julia> C = CalciumField()
Exact Complex Field

julia> exp(C(pi) * C(1im)) + 1
0

julia> log(C(-1))
3.14159*I {a*b where a = 3.14159 [Pi], b = I [b^2+1=0]}

julia> log(C(-1)) ^ 2
-9.86960 {-a^2 where a = 3.14159 [Pi], b = I [b^2+1=0]}

julia> log(C(10)^23) // log(C(100))
11.5000 {23/2}

julia> 4*atan(C(1)//5) - atan(C(1)//239) == C(pi)//4
true

julia> Cx, x = PolynomialRing(C, "x")
(Univariate Polynomial Ring in x over Exact Complex Field, x)

julia> (a, b) = (sqrt(C(2)), sqrt(C(3)))
(1.41421 {a where a = 1.41421 [a^2-2=0]}, 1.73205 {a where a = 1.73205 [a^2-3=0]})

julia> (x-a-b)*(x-a+b)*(x+a-b)*(x+a+b)
x^4 + (-10)*x^2 + 1
```

## Conversions and numerical evaluation

Calcium numbers can created from integers (`ZZ`), rationals (`QQ`)
and algebraic numbers (`QQbar`), and through the application of
arithmetic operations and transcendental functions.

Calcium numbers can be converted to integers, rational and algebraic fields
provided that the values are integer, rational or algebraic.
An exception is thrown if the value does not belong to the target domain,
if Calcium is unable to prove that the value belongs
to the target domain, or if Calcium is unable to compute the explicit
value because of evaluation limits.

```julia
julia> QQ(C(1))
1

julia> QQBar(sqrt(C(2)) // 2)
Root 0.707107 of 2x^2 - 1

julia> QQ(C(pi))
ERROR: unable to convert to a rational number

julia> QQ(C(10) ^ C(10^9))
ERROR: unable to convert to a rational number
```

To compute arbitrary-precision numerical enclosures, convert to
`ArbField` or `AcbField`:

```julia
julia> CC = AcbField(64);

julia> CC(exp(C(1im)))
[0.54030230586813971740 +/- 9.37e-22] + [0.84147098480789650665 +/- 2.51e-21]*im
```

The constructor

```julia
(R::AcbField)(a::ca; parts::Bool=false)
```

returns an enclosure of the complex number *a*.
It attempts to obtain a relative accuracy of *prec* bits where *prec*
is the precision of the target field, but it is not guaranteed that
this goal is achieved.

If *parts* is set to *true*, it attempts to achieve the target
accuracy for both real and imaginary parts. This can be significantly
more expensive if one part is smaller than the other, or if the
number is nontrivially purely real or purely imaginary (in which
case an exact proof attempt is made).

```julia
julia> x = sin(C(1), form=:exponential)
0.841471 + 0e-24*I {(-a^2*b+b)/(2*a) where a = 0.540302 + 0.841471*I [Exp(1.00000*I {b})], b = I [b^2+1=0]}

julia> AcbField(64)(x)
[0.84147098480789650665 +/- 2.51e-21] + [+/- 4.77e-29]*im

julia> AcbField(64)(x, parts=true)
[0.84147098480789650665 +/- 2.51e-21]
```

The constructor

```julia
(R::ArbField)(a::ca; check::Bool=true)
```

returns a real enclosure.
If *check* is set to *true* (default), the number *a* is verified to be
real, and an exception is thrown if this cannot be determined.
With *check* set to *false*, this function returns an enclosure of
the real part of *a* without checking that the imaginary part is zero.
This can be significantly faster.

## Comparisons and properties

Except where otherwise noted, predicate functions such as `iszero`, `==`, `<` and `isreal`
act on the mathematical values of Calcium field elements.
For example, although evaluating $x = \sqrt{2} \sqrt{3}$ and $y = \sqrt{6}$
results in different internal representations ($x \in \mathbb{Q}(\sqrt{3}, \sqrt{2})$
and $y \in \mathbb{Q}(\sqrt{6})$),
the numbers compare as equal:

```julia
julia> x = sqrt(C(2)) * sqrt(C(3))
2.44949 {a*b where a = 1.73205 [a^2-3=0], b = 1.41421 [b^2-2=0]}

julia> y = sqrt(C(6))
2.44949 {a where a = 2.44949 [a^2-6=0]}

julia> x == y
true

julia> iszero(x - y)
true

julia> isinteger(x - y)
true
```

Predicate functions return *true* if the property is provably true
and *false* if the property if provably false. If Calcium is unable
to prove the truth value, an exception is thrown.
For example, with default settings, Calcium is currently
able to prove that $e^{e^{-1000}} \ne 1$,
but it fails to prove $e^{e^{-3000}} \ne 1$:

```julia
julia> x = exp(exp(C(-1000)))
1.00000 {a where a = 1.00000 [Exp(5.07596e-435 {b})], b = 5.07596e-435 [Exp(-1000)]}

julia> x == 1
false

julia> x = exp(exp(C(-3000)))
1.00000 {a where a = 1.00000 [Exp(1.30784e-1303 {b})], b = 1.30784e-1303 [Exp(-3000)]}

julia> x == 1
ERROR: Unable to perform operation (failed deciding truth of a predicate): isequal
...
```

In this case, we can get an answer by allowing a higher working precision:

```julia
julia> C2 = CalciumField(options=Dict(:prec_limit => 10^5));

julia> exp(exp(C2(-3000))) == 1
false
```

Real numbers can be ordered and sorted the usual way.
We illustrate finding square roots that are well-approximated by integers:

```julia
julia> sort([sqrt(C(n)) for n=0:10], by=x -> abs(x - floor(x + C(1)//2)))
11-element Array{ca,1}:
 0
 1
 2
 3
 3.16228 {a where a = 3.16228 [a^2-10=0]}
 2.82843 {2*a where a = 1.41421 [a^2-2=0]}
 2.23607 {a where a = 2.23607 [a^2-5=0]}
 1.73205 {a where a = 1.73205 [a^2-3=0]}
 2.64575 {a where a = 2.64575 [a^2-7=0]}
 1.41421 {a where a = 1.41421 [a^2-2=0]}
 2.44949 {a where a = 2.44949 [a^2-6=0]}
```

As currently implemented, order comparisons involving nonreal numbers yield
*false* (in both directions)
rather than throwing an exception:

```julia
julia> C(1im) < C(1im)
false

julia> C(1im) > C(1im)
false
```

This behavior may be changed or
may become configurable in the future.

**Interface**

```@docs
iszero(a::ca)
isone(a::ca)
isalgebraic(a::ca)
isrational(a::ca)
isinteger(a::ca)
isreal(a::ca)
isimaginary(a::ca)
```

## Infinities and special values

By default, `CalciumField` does not permit creating values that are not
numbers, and any non-number value (unsigned infinity, signed infinity,
Undefined) will result in an exception.
This also applies to the special value Unknown, used in situations
where Calcium is unable to prove that a value is a number.
To enable special values, use `extended=true`.

```julia
julia> C = CalciumField()
Exact Complex Field

julia> 1 // C(0)
ERROR: DomainError with UnsignedInfinity:
Non-number result
...

julia> Cext = CalciumField(extended=true)
Exact Complex Field (Extended)

julia> 1 // Cext(0)
UnsignedInfinity
```

Note that special values do not satisfy the properties of a mathematical
ring or field. You will likely get meaningless results if you put
infinities in matrices or polynomials.

```@docs
isnumber(a::ca)
isundefined(a::ca)
isinf(a::ca)
isuinf(a::ca)
is_signed_inf(a::ca)
isunknown(a::ca)
```

## Complex parts

## Elementary functions

## Special functions

## Rewriting and simplification

```@docs
complex_normal_form(a::ca; deep::Bool=true)
```

