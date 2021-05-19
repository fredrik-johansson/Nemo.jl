```@meta
CurrentModule = Nemo
```

# Algebraic numbers

Nemo allows working with exact real and complex algebraic numbers.

 Library        | Element type  | Parent type
----------------|---------------|--------------------
Calcium         | `qqbar`       | `CalciumQQBarField`

The default algebraic number type in Nemo is provided by Calcium. The
associated field of algebraic numbers is represented by the constant
parent object called `CalciumQQBar`.

The default algebraic number type represents algebraic numbers
in canonical form using minimal polynomials. This works well for representing
individual algebraic numbers, but it does not provide the best
performance for field arithmetic.
To compute in a fixed subfield of $\overline{\mathbb{Q}}$,
it is typically far more efficient to fix a generator $a$
and construct an Antic number field to represent $\mathbb{Q}(a)$.

In the future, Nemo will provide an alternative implementation
of algebraic numbers using Calcium fields, which will construct
efficient internal representations automatically.
The minimal polynomial representation will still be available
since it provides canonical forms with predictable behavior.

