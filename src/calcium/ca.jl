###############################################################################
#
#   ca.jl : Calcium field elements
#
###############################################################################

export ca, CalciumField, inf, uinf, undefined, unknown, const_pi, const_euler,
       const_i, complex_normal_form, csgn,
       gamma, erf, erfi, erfc

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

parent(a::ca) = a.parent

parent_type(::Type{ca}) = CalciumField

elem_type(::Type{CalciumField}) = ca

base_ring(a::CalciumField) = Union{}   #  ?

base_ring(a::ca) = Union{}   #  ?

isdomain_type(::Type{ca}) = true

###############################################################################
#
#   Hashing
#
###############################################################################

# todo: implement nontrivial hash functions on C
function Base.hash(a::ca, h::UInt)
   return h
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::ca) = a

###############################################################################
#
#   I/O
#
###############################################################################

function native_string(x::ca)
   cstr = ccall((:ca_get_str, libcalcium),
        Ptr{UInt8}, (Ref{ca}, Ref{CalciumField}), x, x.parent)
   res = unsafe_string(cstr)
   ccall((:flint_free, libflint), Nothing, (Ptr{UInt8},), cstr)

   return res
end

function show(io::IO, x::ca)
   print(io, native_string(x))
end

needs_parentheses(x::ca) = true

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(C::CalciumField) = C()

function one(C::CalciumField)
   z = ca(C)
   ccall((:ca_one, libcalcium), Nothing,
         (Ref{ca}, Int, Ref{CalciumField}), z, v, C)
   return z
end

###############################################################################
#
#   Random generation
#
###############################################################################

function rand(C::CalciumField; depth::Int, bits::Int,
                                            randtype::Symbol=:null)
   state = _flint_rand_states[Threads.threadid()]
   x = C()

   depth = max(depth, 0)
   bits = max(bits, 1)

   if randtype == :null
      ccall((:ca_randtest, libcalcium), Nothing,
          (Ref{ca}, Ptr{Cvoid}, Int, Int, Ref{CalciumField}),
                x, state.ptr, depth, bits, C)
   elseif randtype == :rational
      ccall((:ca_randtest_rational, libcalcium), Nothing,
          (Ref{ca}, Ptr{Cvoid}, Int, Ref{CalciumField}),
                x, state.ptr, bits, C)
   elseif randtype == :special
      ccall((:ca_randtest, libcalcium), Nothing,
          (Ref{ca}, Ptr{Cvoid}, Int, Int, Ref{CalciumField}),
                x, state.ptr, depth, bits, C)
   else
      error("randtype not defined")
   end

   check_special(x)
   return x
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function _isspecial(a::ca)
   return (a.field & 3) != 0
end

# todo: distiguish unknown
function check_special(a::ca)
   if !a.parent.extended && _isspecial(a)
      throw(DomainError(a, "Non-number result"))
   end
end

function same_parent(a::ca, b::ca)
   if a.parent == b.parent
      return (a, b)
   else
      C = a.parent
      r = C()
      ccall((:ca_transfer, libcalcium), Nothing,
         (Ref{ca}, Ref{CalciumField}, Ref{ca}, Ref{CalciumField}),
         r, a.parent, b, b.parent)
      check_special(r)
      return r
   end
end

function +(a::ca, b::ca)
   a, b = same_parent(a, b)
   C = a.parent
   r = C()
   ccall((:ca_add, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, b, C)
   check_special(r)
   return r
end

function -(a::ca, b::ca)
   a, b = same_parent(a, b)
   C = a.parent
   r = C()
   ccall((:ca_sub, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, b, C)
   check_special(r)
   return r
end

function *(a::ca, b::ca)
   a, b = same_parent(a, b)
   C = a.parent
   r = C()
   ccall((:ca_mul, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, b, C)
   check_special(r)
   return r
end

function //(a::ca, b::ca)
   a, b = same_parent(a, b)
   C = a.parent
   r = C()
   ccall((:ca_div, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, b, C)
   check_special(r)
   return r
end

function ^(a::ca, b::ca)
   a, b = same_parent(a, b)
   C = a.parent
   r = C()
   ccall((:ca_pow, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, b, C)
   check_special(r)
   return r
end

function -(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_neg, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function inv(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_inv, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

###############################################################################
#
#   Comparison and predicates
#
###############################################################################

function ==(a::ca, b::ca)
   a, b = same_parent(a, b)
   C = a.parent
   t = ccall((:ca_check_equal, libcalcium), Cint,
        (Ref{ca}, Ref{ca}, Ref{CalciumField}), a, b, C)
   return truth_as_bool(t, :isequal)
end

function isless(a::ca, b::ca)
   a, b = same_parent(a, b)
   C = a.parent
   t = ccall((:ca_check_lt, libcalcium), Cint,
        (Ref{ca}, Ref{ca}, Ref{CalciumField}), a, b, C)
   return truth_as_bool(t, :isless)
end

isless(a::ca, b::qqbar) = isless(a, parent(a)(b))
isless(a::ca, b::fmpz) = isless(a, parent(a)(b))
isless(a::ca, b::fmpq) = isless(a, parent(a)(b))
isless(a::ca, b::Int) = isless(a, parent(a)(b))
isless(a::qqbar, b::ca) = isless(parent(b)(a), b)
isless(a::fmpq, b::ca) = isless(parent(b)(a), b)
isless(a::fmpz, b::ca) = isless(parent(b)(a), b)
isless(a::Int, b::ca) = isless(parent(b)(a), b)

function isnumber(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_number, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isnumber)
end

function iszero(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_zero, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :iszero)
end

function isone(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_zero, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isone)
end

function isalgebraic(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_algebraic, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isalgebraic)
end

function isrational(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_rational, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isrational)
end

function isinteger(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_integer, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isinteger)
end

function isreal(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_real, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isreal)
end

function isimaginary(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_imaginary, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isimaginary)
end

function isundefined(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_undefined, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isundefined)
end

function isinf(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_infinity, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isinf)
end

function isuinf(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_uinf, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :isuinf)
end

function is_signed_inf(a::ca)
   C = a.parent
   t = ccall((:ca_check_is_signed_inf, libcalcium), Cint,
        (Ref{ca}, Ref{CalciumField}), a, C)
   return truth_as_bool(t, :is_signed_inf)
end

###############################################################################
#
#   Special values and constants
#
###############################################################################

function const_pi(C::CalciumField)
   r = C()
   ccall((:ca_pi, libcalcium), Nothing, (Ref{ca}, Ref{CalciumField}), r, C)
   return r
end

function const_euler(C::CalciumField)
   r = C()
   ccall((:ca_euler, libcalcium), Nothing, (Ref{ca}, Ref{CalciumField}), r, C)
   return r
end

function const_i(C::CalciumField)
   r = C()
   ccall((:ca_i, libcalcium), Nothing, (Ref{ca}, Ref{CalciumField}), r, C)
   return r
end

function uinf(C::CalciumField)
   r = C()
   ccall((:ca_uinf, libcalcium), Nothing,
         (Ref{ca}, Ref{CalciumField}), r, C)
   check_special(r)
   return r
end

function inf(C::CalciumField)
   r = C()
   ccall((:ca_pos_inf, libcalcium), Nothing,
         (Ref{ca}, Ref{CalciumField}), r, C)
   check_special(r)
   return r
end

function undefined(C::CalciumField)
   r = C()
   ccall((:ca_undefined, libcalcium), Nothing,
         (Ref{ca}, Ref{CalciumField}), r, C)
   check_special(r)
   return r
end

function unknown(C::CalciumField)
   r = C()
   ccall((:ca_unknown, libcalcium), Nothing,
         (Ref{ca}, Ref{CalciumField}), r, C)
   check_special(r)
   return r
end

###############################################################################
#
#   Complex parts
#
###############################################################################

function real(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_re, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function imag(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_im, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function angle(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_arg, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function csgn(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_csgn, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function sign(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_sgn, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function abs(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_abs, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function conj(a::ca; form::Symbol=:default)
   C = a.parent
   r = C()
   if form == :default
      ccall((:ca_conj, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :deep
      ccall((:ca_conj_deep, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :shallow
      ccall((:ca_conj_shallow, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   else
      error("unknown form: ", form)
   end
   check_special(r)
   return r
end

function floor(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_floor, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function ceil(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_ceil, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

###############################################################################
#
#   Elementary functions
#
###############################################################################

function sqrt(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_sqrt, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function exp(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_exp, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function log(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_log, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function sin(a::ca; form::Symbol=:default)
   C = a.parent
   r = C()
   if form == :default
      ccall((:ca_sin, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :exponential
      ccall((:ca_sin_cos_exponential, libcalcium), Nothing,
             (Ref{ca}, Ptr{Nothing}, Ref{ca}, Ref{CalciumField}), r, C_NULL, a, C)
   elseif form == :tangent
      ccall((:ca_sin_cos_tangent, libcalcium), Nothing,
             (Ref{ca}, Ptr{Nothing}, Ref{ca}, Ref{CalciumField}), r, C_NULL, a, C)
   elseif form == :direct
      ccall((:ca_sin_cos_direct, libcalcium), Nothing,
             (Ref{ca}, Ptr{Nothing}, Ref{ca}, Ref{CalciumField}), r, C_NULL, a, C)
   else
      error("unknown form: ", form)
   end
   check_special(r)
   return r
end

function cos(a::ca; form::Symbol=:default)
   C = a.parent
   r = C()
   if form == :default
      ccall((:ca_cos, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :exponential
      ccall((:ca_sin_cos_exponential, libcalcium), Nothing,
             (Ptr{Nothing}, Ref{ca}, Ref{ca}, Ref{CalciumField}), C_NULL, r, a, C)
   elseif form == :tangent
      ccall((:ca_sin_cos_tangent, libcalcium), Nothing,
             (Ptr{Nothing}, Ref{ca}, Ref{ca}, Ref{CalciumField}), C_NULL, r, a, C)
   elseif form == :direct || form == :sine_cosine
      ccall((:ca_sin_cos_direct, libcalcium), Nothing,
             (Ptr{Nothing}, Ref{ca}, Ref{ca}, Ref{CalciumField}), C_NULL, r, a, C)
   else
      error("unknown form: ", form)
   end
   check_special(r)
   return r
end

function tan(a::ca; form::Symbol=:default)
   C = a.parent
   r = C()
   if form == :default
      ccall((:ca_tan, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :exponential
      ccall((:ca_tan_exponential, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :direct || form == :tangent
      ccall((:ca_tan_direct, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :sine_cosine
      ccall((:ca_tan_sine_cosine, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   else
      error("unknown form: ", form)
   end
   check_special(r)
   return r
end

function atan(a::ca; form::Symbol=:default)
   C = a.parent
   r = C()
   if form == :default
      ccall((:ca_atan, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :logarithm
      ccall((:ca_atan_logarithm, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :direct || form == :arctangent
      ccall((:ca_atan_direct, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   else
      error("unknown form: ", form)
   end
   check_special(r)
   return r
end

function asin(a::ca; form::Symbol=:default)
   C = a.parent
   r = C()
   if form == :default
      ccall((:ca_asin, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :logarithm
      ccall((:ca_asin_logarithm, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :direct
      ccall((:ca_asin_direct, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   else
      error("unknown form: ", form)
   end
   check_special(r)
   return r
end

function acos(a::ca; form::Symbol=:default)
   C = a.parent
   r = C()
   if form == :default
      ccall((:ca_acos, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :logarithm
      ccall((:ca_acos_logarithm, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   elseif form == :direct
      ccall((:ca_acos_direct, libcalcium), Nothing,
             (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   else
      error("unknown form: ", form)
   end
   check_special(r)
   return r
end

function gamma(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_gamma, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function erf(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_erf, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function erfi(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_erfi, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

function erfc(a::ca)
   C = a.parent
   r = C()
   ccall((:ca_erfc, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{CalciumField}), r, a, C)
   check_special(r)
   return r
end

###############################################################################
#
#   Rewriting and normal forms
#
###############################################################################

@doc Markdown.doc"""
    complex_normal_form(a::ca, deep::Bool=true)

Returns the input rewritten using standardizing transformations over the
complex numbers:

* Elementary functions are rewritten in terms of exponentials, roots
  and logarithms.

* Complex parts are rewritten using logarithms, square roots, and (deep)
  complex conjugates.

* Algebraic numbers are rewritten in terms of cyclotomic fields where
  applicable.

If deep is set, the rewriting is applied recursively to the tower of
extension numbers; otherwise, the rewriting is only applied to the
top-level extension numbers.

The result is not a normal form in the strong sense (the same number
can have many possible representations even after applying this
transformation), but this transformation can nevertheless be a useful
heuristic for simplification.
"""
function complex_normal_form(a::ca; deep::Bool=true)
   C = a.parent
   r = C()
   ccall((:ca_rewrite_complex_normal_form, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Cint, Ref{CalciumField}), r, a, deep, C)
   check_special(r)
   return r
end

###############################################################################
#
#   Conversions
#
###############################################################################

# todo: wrap qqbar_get_fmpq, qqbar_get_fmpz
# todo: field conversion versions

function fmpq(a::ca)
   C = a.parent
   res = fmpq()
   ok = Bool(ccall((:ca_get_fmpq, libcalcium), Cint,
        (Ref{fmpq}, Ref{ca}, Ref{CalciumField}), res, a, C))
   !ok && error("unable to convert to a rational number")
   return res
end

function fmpz(a::ca)
   C = a.parent
   res = fmpz()
   ok = Bool(ccall((:ca_get_fmpz, libcalcium), Cint,
        (Ref{fmpz}, Ref{ca}, Ref{CalciumField}), res, a, C))
   !ok && error("unable to convert to an integer")
   return res
end

function qqbar(a::ca)
   C = a.parent
   res = qqbar()
   ok = Bool(ccall((:ca_get_qqbar, libcalcium), Cint,
        (Ref{qqbar}, Ref{ca}, Ref{CalciumField}), res, a, C))
   !ok && error("unable to convert to an algebraic number")
   return res
end

(R::FlintRationalField)(a::ca) = fmpq(a)
(R::FlintIntegerRing)(a::ca) = fmpz(a)
(R::CalciumQQBarField)(a::ca) = qqbar(a)

function (R::AcbField)(a::ca; parts::Bool=false)
   C = a.parent
   prec = precision(R)
   z = R()
   if parts
      ccall((:ca_get_acb_accurate_parts, libcalcium),
        Nothing, (Ref{acb}, Ref{ca}, Int, Ref{CalciumField}), z, a, prec, C)
   else
      ccall((:ca_get_acb, libcalcium),
        Nothing, (Ref{acb}, Ref{ca}, Int, Ref{CalciumField}), z, a, prec, C)
   end
   return z
end

function (R::ArbField)(a::ca; check::Bool=true)
   C = a.parent
   prec = precision(R)
   if check
      z = AcbField(prec)(a, parts=true)
      if isreal(z)
         return real(z)
      else
         error("unable to convert to a real number")
      end
   else
      z = AcbField(prec)(a, parts=false)
      if accuracy_bits(z) < prec - 5
          z = AcbField(prec)(a, parts=true)
      end
      return real(z)
   end
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(z::ca)
   C = z.parent
   ccall((:ca_zero, libcalcium), Nothing, (Ref{ca}, Ref{CalciumField}), z, C)
   return z
end

function mul!(z::ca, x::ca, y::ca)
   if z.parent != x.parent || x.parent != y.parent
      error("different parents in in-place operation")
   end
   C = z.parent
   ccall((:ca_mul, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{ca}, Ref{CalciumField}), z, x, y, C)
   check_special(z)
   return z
end

function addeq!(z::ca, x::ca)
   if z.parent != x.parent
      error("different parents in in-place operation")
   end
   C = z.parent
   ccall((:ca_add, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{ca}, Ref{CalciumField}), z, z, x, C)
   check_special(z)
   return z
end

function add!(z::ca, x::ca, y::ca)
   if z.parent != x.parent || x.parent != y.parent
      error("different parents in in-place operation")
   end
   C = z.parent
   ccall((:ca_add, libcalcium), Nothing,
         (Ref{ca}, Ref{ca}, Ref{ca}, Ref{CalciumField}), z, x, y, C)
   check_special(z)
   return z
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (C::CalciumField)()
   z = ca(C)
   return z
end

function (C::CalciumField)(v::ca)
   D = v.parent
   if C == D
      return v
   end
   r = C()
   ccall((:ca_transfer, libcalcium), Nothing,
      (Ref{ca}, Ref{CalciumField}, Ref{ca}, Ref{CalciumField}),
      r, C, v, D)
   check_special(r)
   return r
end

function (C::CalciumField)(v::Int)
   z = ca(C)
   ccall((:ca_set_si, libcalcium), Nothing,
         (Ref{ca}, Int, Ref{CalciumField}), z, v, C)
   return z
end

function (C::CalciumField)(v::fmpz)
   z = ca(C)
   ccall((:ca_set_fmpz, libcalcium), Nothing,
         (Ref{ca}, Ref{fmpz}, Ref{CalciumField}), z, v, C)
   return z
end

function (C::CalciumField)(v::fmpq)
   z = ca(C)
   ccall((:ca_set_fmpq, libcalcium), Nothing,
         (Ref{ca}, Ref{fmpq}, Ref{CalciumField}), z, v, C)
   return z
end

function (C::CalciumField)(v::qqbar)
   z = ca(C)
   ccall((:ca_set_qqbar, libcalcium), Nothing,
         (Ref{ca}, Ref{qqbar}, Ref{CalciumField}), z, v, C)
   return z
end

# todo: optimize
function (C::CalciumField)(v::Complex{Int})
   return C(QQBar(v))
end

function (C::CalciumField)(x::Irrational)
  if x == pi
    return const_pi(C)
  else
    error("constant not supported")
  end
end

