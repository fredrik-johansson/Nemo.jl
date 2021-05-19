###############################################################################
#
#   qqbar.jl : Calcium algebraic numbers in minimal polynomial representation
#
###############################################################################

export qqbar, CalciumQQBar, CalciumQQBarField, is_algebraic_integer

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

parent(a::qqbar) = CalciumQQBar

parent_type(::Type{qqbar}) = CalciumQQBarField

elem_type(::Type{CalciumQQBarField}) = qqbar

base_ring(a::CalciumQQBarField) = CalciumQQBar

base_ring(a::qqbar) = CalciumQQBar

isdomain_type(::Type{qqbar}) = true

###############################################################################
#
#   Hashing
#
###############################################################################

# todo
# function Base.hash(a::qqbar, h::UInt)
# end

###############################################################################
#
#   Constructors
#
###############################################################################

function qqbar(a::Int)
   z = qqbar()
   ccall((:qqbar_set_si, libcalcium), Nothing, (Ref{qqbar}, Int, ), z, a)
  return z
end

function qqbar(a::fmpz)
   z = qqbar()
   ccall((:qqbar_set_fmpz, libcalcium), Nothing, (Ref{qqbar}, Ref{fmpz}, ), z, a)
   return z
end

function qqbar(a::fmpq)
   z = qqbar()
   ccall((:qqbar_set_fmpq, libcalcium), Nothing, (Ref{qqbar}, Ref{fmpq}, ), z, a)
   return z
end


###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::qqbar) = a

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

# todo
# function expressify(a::qqbar; context = nothing)::Any
# end

function native_string(x::qqbar)
   cstr = ccall((:qqbar_get_str_nd, libcalcium), Ptr{UInt8},
                (Ref{qqbar}, Int), x, Int(6))
   res = unsafe_string(cstr)
   ccall((:flint_free, libflint), Nothing,
         (Ptr{UInt8},),
         cstr)
   return res
end

function show(io::IO, F::CalciumQQBarField)
  print(io, "Field of algebraic numbers in minimal polynomial representation")
end

function show(io::IO, x::qqbar)
   print(io, native_string(x))
end

needs_parentheses(x::qqbar) = false

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function degree(x::qqbar)
   return ccall((:qqbar_degree, libcalcium), Int, (Ref{qqbar}, ), x)
end

function iszero(x::qqbar)
   return Bool(ccall((:qqbar_is_zero, libcalcium), Cint, (Ref{qqbar},), x))
end

function isone(x::qqbar)
   return Bool(ccall((:qqbar_is_one, libcalcium), Cint, (Ref{qqbar},), x))
end

function isinteger(x::qqbar)
   return Bool(ccall((:qqbar_is_integer, libcalcium), Cint, (Ref{qqbar},), x))
end

function isrational(x::qqbar)
   return Bool(ccall((:qqbar_is_rational, libcalcium), Cint, (Ref{qqbar},), x))
end

function isreal(x::qqbar)
   return Bool(ccall((:qqbar_is_real, libcalcium), Cint, (Ref{qqbar},), x))
end

function is_algebraic_integer(x::qqbar)
   return Bool(ccall((:qqbar_is_algebraic_integer, libcalcium), Cint, (Ref{qqbar},), x))
end

function minpoly(R::FmpzPolyRing, x::qqbar)
   z = R()
   ccall((:fmpz_poly_set, libflint), Nothing, (Ref{fmpz_poly}, Ref{qqbar}, ), z, x)
   return z
end

function minpoly(R::FmpqPolyRing, x::qqbar)
   z = R()
   ccall((:fmpq_poly_set_fmpz_poly, libflint), Nothing, (Ref{fmpq_poly}, Ref{qqbar}, ), z, x)
   return z
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::qqbar)
   z = qqbar()
   ccall((:qqbar_neg, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}), z, a)
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::qqbar, b::qqbar)
   z = qqbar()
   ccall((:qqbar_add, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, a, b)
   return z
end

function -(a::qqbar, b::qqbar)
   z = qqbar()
   ccall((:qqbar_sub, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, a, b)
   return z
end

function *(a::qqbar, b::qqbar)
   z = qqbar()
   ccall((:qqbar_mul, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, a, b)
   return z
end

function ^(a::qqbar, b::qqbar)
   z = qqbar()
   ok = Bool(ccall((:qqbar_pow, libcalcium), Cint,
         (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, a, b))
   !ok && throw(DomainError)
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::qqbar, b::qqbar)
   iszero(b) && throw(DivideError())
   z = qqbar()
   ccall((:qqbar_div, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, a, b)
   return z
end

div(a::qqbar, b::qqbar) = divexact(a, b)

###############################################################################
#
#   Roots
#
###############################################################################

function sqrt(a::qqbar)
   z = qqbar()
   ccall((:qqbar_sqrt, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}), z, a)
   return z
end

function root(a::qqbar, n::Int)
   n <= 0 && throw(DomainError())
   z = qqbar()
   ccall((:qqbar_root_ui, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}, UInt), z, a, n)
   return z
end

function root_of_unity(C::CalciumQQBarField, n::Int)
   n <= 0 && throw(DomainError())
   z = qqbar()
   ccall((:qqbar_root_of_unity, libcalcium), Nothing, (Ref{qqbar}, Int, UInt), z, 1, n)
   return z
end

function root_of_unity(C::CalciumQQBarField, n::Int, k::Int)
   n <= 0 && throw(DomainError())
   z = qqbar()
   ccall((:qqbar_root_of_unity, libcalcium), Nothing, (Ref{qqbar}, Int, UInt), z, k, n)
   return z
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

(a::CalciumQQBarField)() = qqbar()

(a::CalciumQQBarField)(b::Int) = qqbar(b)

(a::CalciumQQBarField)(b::fmpz) = qqbar(b)

(a::CalciumQQBarField)(b::fmpq) = qqbar(b)

