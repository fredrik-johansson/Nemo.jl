```@meta
CurrentModule = Nemo
```

# Exact real and complex numbers

Exact real and complex numbers are provided by Calcium.
Internally, a number $z$ is represented as an element
of an extension field of the rational numbers. That is,

$z \in \mathbb{Q}(a_1,\ldots,a_n)$

where $a_1, \ldots, a_n$
are symbolically defined
real or complex numbers such as $\pi$, $\sqrt{2}$ or $e^{\sqrt{2} \pi i}$.
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


