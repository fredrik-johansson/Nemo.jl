###############################################################################
#
#   ca.jl : Calcium field elements
#
###############################################################################

export ca, CalciumField, inf, uinf, undefined, unknown, const_pi, const_euler,
       const_i

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

###############################################################################
#
#   Constructors
#
###############################################################################

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

###############################################################################
#
#   Random generation
#
###############################################################################

###############################################################################
#
#   Arithmetic
#
###############################################################################

function _isspecial(a::ca)
   return (a.field & 3) != 0
end

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

###############################################################################
#
#   Comparison
#
###############################################################################

###############################################################################
#
#   Complex parts
#
###############################################################################

###############################################################################
#
#   Elementary functions
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

#=
function sin(a::ca, form::Symbol=:default)
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
=#

###############################################################################
#
#   Conversions
#
###############################################################################

###############################################################################
#
#   Unsafe functions
#
###############################################################################

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (C::CalciumField)()
   z = ca(C)
   return z
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

