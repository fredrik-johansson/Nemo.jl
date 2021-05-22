###############################################################################
#
#   qqbar.jl : Calcium algebraic numbers in minimal polynomial representation
#
###############################################################################

export qqbar, CalciumQQBar, CalciumQQBarField, is_algebraic_integer, rand, abs2,
       csgn, sign_real, sign_imag, fmpq, fmpz, exp_pi_i, atanpi, asinpi, acospi,
       conjugates, eigenvalues, guess, root_of_unity_as_args, is_root_of_unity,
       log_pi_i

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

# todo: want a C function for this
function Base.hash(a::qqbar, h::UInt)
   R, x = PolynomialRing(FlintZZ, "x")
   return xor(hash(minpoly(R, a)), h)
end

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
  print(io, "Field of Algebraic Numbers in minimal polynomial representation")
end

function show(io::IO, x::qqbar)
   print(io, native_string(x))
end

needs_parentheses(x::qqbar) = true

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(a::CalciumQQBarField) = a(0)

one(a::CalciumQQBarField) = a(1)

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
#   Random generation
#
###############################################################################

function rand(R::CalciumQQBarField; degree::Int, bits::Int, randtype::Symbol=:null)
   state = _flint_rand_states[Threads.threadid()]
   x = R()

   degree <= 0 && error("degree must be positive")
   bits <= 0 && error("bits must be positive")

   if randtype == :null
      ccall((:qqbar_randtest, libcalcium), Nothing,
          (Ref{qqbar}, Ptr{Cvoid}, Int, Int), x, state.ptr, degree, bits)
   elseif randtype == :real
      ccall((:qqbar_randtest_real, libcalcium), Nothing,
          (Ref{qqbar}, Ptr{Cvoid}, Int, Int), x, state.ptr, degree, bits)
   elseif randtype == :nonreal
      degree < 2 && error("nonreal requires degree >= 2")
      ccall((:qqbar_randtest_nonreal, libcalcium), Nothing,
          (Ref{qqbar}, Ptr{Cvoid}, Int, Int), x, state.ptr, degree, bits)
   else
      error("randtype not defined")
   end

   return x
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

function +(a::qqbar, b::fmpq)
   z = qqbar()
   ccall((:qqbar_add_fmpq, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{fmpq}), z, a, b)
   return z
end

function +(a::qqbar, b::fmpz)
   z = qqbar()
   ccall((:qqbar_add_fmpz, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{fmpz}), z, a, b)
   return z
end

function +(a::qqbar, b::Int)
   z = qqbar()
   ccall((:qqbar_add_si, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Int), z, a, b)
   return z
end

+(a::fmpq, b::qqbar) = b + a
+(a::fmpz, b::qqbar) = b + a
+(a::Int, b::qqbar) = b + a

function -(a::qqbar, b::qqbar)
   z = qqbar()
   ccall((:qqbar_sub, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, a, b)
   return z
end

function -(a::qqbar, b::fmpq)
   z = qqbar()
   ccall((:qqbar_sub_fmpq, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{fmpq}), z, a, b)
   return z
end

function -(a::qqbar, b::fmpz)
   z = qqbar()
   ccall((:qqbar_sub_fmpz, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{fmpz}), z, a, b)
   return z
end

function -(a::qqbar, b::Int)
   z = qqbar()
   ccall((:qqbar_sub_si, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Int), z, a, b)
   return z
end

function -(a::fmpq, b::qqbar)
   z = qqbar()
   ccall((:qqbar_fmpq_sub, libcalcium), Nothing,
         (Ref{qqbar}, Ref{fmpq}, Ref{qqbar}), z, a, b)
   return z
end

function -(a::fmpz, b::qqbar)
   z = qqbar()
   ccall((:qqbar_fmpz_sub, libcalcium), Nothing,
         (Ref{qqbar}, Ref{fmpz}, Ref{qqbar}), z, a, b)
   return z
end

function -(a::Int, b::qqbar)
   z = qqbar()
   ccall((:qqbar_si_sub, libcalcium), Nothing,
         (Ref{qqbar}, Int, Ref{qqbar}), z, a, b)
   return z
end

function *(a::qqbar, b::qqbar)
   z = qqbar()
   ccall((:qqbar_mul, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, a, b)
   return z
end

function *(a::qqbar, b::fmpq)
   z = qqbar()
   ccall((:qqbar_mul_fmpq, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{fmpq}), z, a, b)
   return z
end

function *(a::qqbar, b::fmpz)
   z = qqbar()
   ccall((:qqbar_mul_fmpz, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{fmpz}), z, a, b)
   return z
end

function *(a::qqbar, b::Int)
   z = qqbar()
   ccall((:qqbar_mul_si, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Int), z, a, b)
   return z
end

*(a::fmpq, b::qqbar) = b * a
*(a::fmpz, b::qqbar) = b * a
*(a::Int, b::qqbar) = b * a

function ^(a::qqbar, b::qqbar)
   z = qqbar()
   ok = Bool(ccall((:qqbar_pow, libcalcium), Cint,
         (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, a, b))
   !ok && throw(DomainError((a, b)))
   return z
end

# todo: want qqbar_pow_fmpz, qqbar_pow_fmpq, qqbar_pow_si
^(a::qqbar, b::fmpz) = a ^ qqbar(b)
^(a::qqbar, b::fmpq) = a ^ qqbar(b)
^(a::qqbar, b::Int) = a ^ qqbar(b)
^(a::fmpz, b::qqbar) = qqbar(a) ^ b
^(a::fmpq, b::qqbar) = qqbar(a) ^ b
^(a::Int, b::qqbar) = qqbar(a) ^ b

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

function divexact(a::qqbar, b::fmpq)
   iszero(b) && throw(DivideError())
   z = qqbar()
   ccall((:qqbar_div_fmpq, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{fmpq}), z, a, b)
   return z
end

function divexact(a::qqbar, b::fmpz)
   iszero(b) && throw(DivideError())
   z = qqbar()
   ccall((:qqbar_div_fmpz, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Ref{fmpz}), z, a, b)
   return z
end

function divexact(a::qqbar, b::Int)
   iszero(b) && throw(DivideError())
   z = qqbar()
   ccall((:qqbar_div_si, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Int), z, a, b)
   return z
end

//(a::qqbar, b::qqbar) = divexact(a, b)
//(a::qqbar, b::fmpq) = divexact(a, b)
//(a::qqbar, b::fmpz) = divexact(a, b)
//(a::qqbar, b::Int) = divexact(a, b)
//(a::fmpq, b::qqbar) = divexact(a, b)
//(a::fmpz, b::qqbar) = divexact(a, b)
//(a::Int, b::qqbar) = divexact(a, b)


function <<(a::qqbar, b::Int)
   z = qqbar()
   ccall((:qqbar_mul_2exp_si, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Int), z, a, b)
   return z
end

function >>(a::qqbar, b::Int)
   z = qqbar()
   ccall((:qqbar_mul_2exp_si, libcalcium), Nothing,
         (Ref{qqbar}, Ref{qqbar}, Int), z, a, -b)
   return z
end


###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::qqbar, b::qqbar)
   return Bool(ccall((:qqbar_equal, libcalcium), Cint,
                (Ref{qqbar}, Ref{qqbar}), a, b))
end

function cmp(a::qqbar, b::qqbar)
   !isreal(a) && throw(DomainError(a, "comparing nonreal numbers"))
   !isreal(b) && throw(DomainError(b, "comparing nonreal numbers"))
   return ccall((:qqbar_cmp_re, libcalcium), Cint,
                (Ref{qqbar}, Ref{qqbar}), a, b)
end

function isless(a::qqbar, b::qqbar)
    return cmp(a, b) < 0
end

isless(a::qqbar, b::fmpz) = isless(a, qqbar(b))
isless(a::qqbar, b::fmpq) = isless(a, qqbar(b))
isless(a::qqbar, b::Int) = isless(a, qqbar(b))
isless(a::fmpq, b::qqbar) = isless(qqbar(a), b)
isless(a::fmpz, b::qqbar) = isless(qqbar(a), b)
isless(a::Int, b::qqbar) = isless(qqbar(a), b)

# todo: name and export the following functions?

function cmp_real(a::qqbar, b::qqbar)
   return ccall((:qqbar_cmp_re, libcalcium), Cint,
                (Ref{qqbar}, Ref{qqbar}), a, b)
end

function cmp_imag(a::qqbar, b::qqbar)
   return ccall((:qqbar_cmp_im, libcalcium), Cint,
                (Ref{qqbar}, Ref{qqbar}), a, b)
end

function cmpabs_real(a::qqbar, b::qqbar)
   return ccall((:qqbar_cmpabs_re, libcalcium), Cint,
                (Ref{qqbar}, Ref{qqbar}), a, b)
end

function cmpabs_imag(a::qqbar, b::qqbar)
   return ccall((:qqbar_cmpabs_im, libcalcium), Cint,
                (Ref{qqbar}, Ref{qqbar}), a, b)
end

function cmp_root_order(a::qqbar, b::qqbar)
   return ccall((:qqbar_cmp_root_order, libcalcium), Cint,
                (Ref{qqbar}, Ref{qqbar}), a, b)
end

function isless_root_order(a::qqbar, b::qqbar)
    return cmp_root_order(a, b) < 0
end

# todo: wrap qqbar_equal_fmpq_poly_val

###############################################################################
#
#   Complex parts
#
###############################################################################

function real(a::qqbar)
   z = qqbar()
   ccall((:qqbar_re, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}), z, a)
   return z
end

function imag(a::qqbar)
   z = qqbar()
   ccall((:qqbar_im, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}), z, a)
   return z
end

function abs(a::qqbar)
   z = qqbar()
   ccall((:qqbar_abs, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}), z, a)
   return z
end

function conj(a::qqbar)
   z = qqbar()
   ccall((:qqbar_conj, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}), z, a)
   return z
end

function abs2(a::qqbar)
   z = qqbar()
   ccall((:qqbar_abs2, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}), z, a)
   return z
end

function sign(a::qqbar)
   z = qqbar()
   ccall((:qqbar_sgn, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}), z, a)
   return z
end

function csgn(a::qqbar)
   return qqbar(Int(ccall((:qqbar_csgn, libcalcium), Cint, (Ref{qqbar}, ), a)))
end

function sign_real(a::qqbar)
   return qqbar(Int(ccall((:qqbar_sgn_re, libcalcium), Cint, (Ref{qqbar}, ), a)))
end

function sign_imag(a::qqbar)
   return qqbar(Int(ccall((:qqbar_sgn_im, libcalcium), Cint, (Ref{qqbar}, ), a)))
end

function floor(a::qqbar)
   z = fmpz()
   ccall((:qqbar_floor, libcalcium), Nothing, (Ref{fmpz}, Ref{qqbar}, ), z, a)
   return qqbar(z)
end

function ceil(a::qqbar)
   z = fmpz()
   ccall((:qqbar_ceil, libcalcium), Nothing, (Ref{fmpz}, Ref{qqbar}, ), z, a)
   return qqbar(z)
end


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
   n <= 0 && throw(DomainError(n))
   z = qqbar()
   ccall((:qqbar_root_ui, libcalcium), Nothing, (Ref{qqbar}, Ref{qqbar}, UInt), z, a, n)
   return z
end

function qqbar_vec(n::Int)
   return ccall((:_qqbar_vec_init, libcalcium), Ptr{qqbar_struct}, (Int,), n)
end

function array(R::CalciumQQBarField, v::Ptr{qqbar_struct}, n::Int)
   r = Vector{qqbar}(undef, n)
   for i=1:n
       r[i] = R()
       ccall((:qqbar_set, libcalcium), Nothing, (Ref{qqbar}, Ptr{qqbar_struct}),
           r[i], v + (i-1)*sizeof(qqbar_struct))
   end
   return r
end

function qqbar_vec_clear(v::Ptr{qqbar_struct}, n::Int)
   ccall((:_qqbar_vec_clear, libcalcium), Nothing, (Ptr{qqbar_struct}, Int), v, n)
end

function roots(f::fmpz_poly, R::CalciumQQBarField)
   deg = degree(f)
   if deg <= 0
      return Array{qqbar}(undef, 0)
   end
   roots = qqbar_vec(deg)
   ccall((:qqbar_roots_fmpz_poly, libcalcium), Nothing, (Ptr{qqbar_struct}, Ref{fmpz_poly}, Int), roots, f, 0)
   res = array(R, roots, deg)
   qqbar_vec_clear(roots, deg)
   return res
end

function roots(f::fmpq_poly, R::CalciumQQBarField)
   deg = degree(f)
   if deg <= 0
      return Array{qqbar}(undef, 0)
   end
   roots = qqbar_vec(deg)
   ccall((:qqbar_roots_fmpq_poly, libcalcium), Nothing, (Ptr{qqbar_struct}, Ref{fmpq_poly}, Int), roots, f, 0)
   res = array(R, roots, deg)
   qqbar_vec_clear(roots, deg)
   return res
end

function conjugates(a::qqbar)
   deg = degree(a)
   if deg == 1
      return [a]
   end
   conjugates = qqbar_vec(deg)
   ccall((:qqbar_conjugates, libcalcium), Nothing, (Ptr{qqbar_struct}, Ref{qqbar}), conjugates, a)
   res = array(parent(a), conjugates, deg)
   qqbar_vec_clear(conjugates, deg)
   return res
end

function eigenvalues(A::fmpz_mat, R::CalciumQQBarField)
   n = nrows(A)
   !issquare(A) && throw(DomainError(A, "a square matrix is required"))
   if n == 0
      return Array{qqbar}(undef, 0)
   end
   roots = qqbar_vec(n)
   ccall((:qqbar_eigenvalues_fmpz_mat, libcalcium), Nothing, (Ptr{qqbar_struct}, Ref{fmpz_mat}, Int), roots, A, 0)
   res = array(R, roots, n)
   qqbar_vec_clear(roots, n)
   return res
end

function eigenvalues(A::fmpq_mat, R::CalciumQQBarField)
   n = nrows(A)
   !issquare(A) && throw(DomainError(A, "a square matrix is required"))
   if n == 0
      return Array{qqbar}(undef, 0)
   end
   roots = qqbar_vec(n)
   ccall((:qqbar_eigenvalues_fmpq_mat, libcalcium), Nothing, (Ptr{qqbar_struct}, Ref{fmpq_mat}, Int), roots, A, 0)
   res = array(R, roots, n)
   qqbar_vec_clear(roots, n)
   return res
end

###############################################################################
#
#   Roots of unity and trigonometric functions
#
###############################################################################

function root_of_unity(C::CalciumQQBarField, n::Int)
   n <= 0 && throw(DomainError(n))
   z = qqbar()
   ccall((:qqbar_root_of_unity, libcalcium), Nothing, (Ref{qqbar}, Int, UInt), z, 1, n)
   return z
end

function root_of_unity(C::CalciumQQBarField, n::Int, k::Int)
   n <= 0 && throw(DomainError(n))
   z = qqbar()
   ccall((:qqbar_root_of_unity, libcalcium), Nothing, (Ref{qqbar}, Int, UInt), z, k, n)
   return z
end

function is_root_of_unity(a::qqbar)
   return Bool(ccall((:qqbar_is_root_of_unity, libcalcium), Cint, (Ptr{Int}, Ptr{Int}, Ref{qqbar}), C_NULL, C_NULL, a))
end

function root_of_unity_as_args(a::qqbar)
   p = Vector{Int}(undef, 1)
   q = Vector{Int}(undef, 1)
   if !Bool(ccall((:qqbar_is_root_of_unity, libcalcium), Cint, (Ptr{Int}, Ptr{Int}, Ref{qqbar}), p, q, a))
      throw(DomainError(a, "value is not a root of unity"))
   end
   return (q[1], p[1])
end

function exp_pi_i(a::qqbar)
   r = fmpq(a)
   p = Int(numerator(r))
   q = Int(denominator(r))
   z = qqbar()
   ccall((:qqbar_exp_pi_i, libcalcium), Nothing, (Ref{qqbar}, Int, Int), z, p, q)
   return z
end

function sinpi(a::qqbar)
   r = fmpq(a)
   p = Int(numerator(r))
   q = Int(denominator(r))
   z = qqbar()
   ccall((:qqbar_sin_pi, libcalcium), Nothing, (Ref{qqbar}, Int, Int), z, p, q)
   return z
end

function cospi(a::qqbar)
   r = fmpq(a)
   p = Int(numerator(r))
   q = Int(denominator(r))
   z = qqbar()
   ccall((:qqbar_cos_pi, libcalcium), Nothing, (Ref{qqbar}, Int, Int), z, p, q)
   return z
end

function tanpi(a::qqbar)
   r = fmpq(a)
   p = Int(numerator(r))
   q = Int(denominator(r))
   z = qqbar()
   if !Bool(ccall((:qqbar_tan_pi, libcalcium), Cint, (Ref{qqbar}, Int, Int), z, p, q))
      throw(DomainError(a, "function value is not algebraic"))
   end
   return z
end

function atanpi(a::qqbar)
   p = Vector{Int}(undef, 1)
   q = Vector{Int}(undef, 1)
   if !Bool(ccall((:qqbar_atan_pi, libcalcium), Cint, (Ptr{Int}, Ptr{Int}, Ref{qqbar}), p, q, a))
      throw(DomainError(a, "function value is not algebraic"))
   end
   return qqbar(p[1]) // q[1]
end

function asinpi(a::qqbar)
   p = Vector{Int}(undef, 1)
   q = Vector{Int}(undef, 1)
   if !Bool(ccall((:qqbar_asin_pi, libcalcium), Cint, (Ptr{Int}, Ptr{Int}, Ref{qqbar}), p, q, a))
      throw(DomainError(a, "function value is not algebraic"))
   end
   return qqbar(p[1]) // q[1]
end

function acospi(a::qqbar)
   p = Vector{Int}(undef, 1)
   q = Vector{Int}(undef, 1)
   if !Bool(ccall((:qqbar_acos_pi, libcalcium), Cint, (Ptr{Int}, Ptr{Int}, Ref{qqbar}), p, q, a))
      throw(DomainError(a, "function value is not algebraic"))
   end
   return qqbar(p[1]) // q[1]
end

function log_pi_i(a::qqbar)
   p = Vector{Int}(undef, 1)
   q = Vector{Int}(undef, 1)
   if !Bool(ccall((:qqbar_log_pi_i, libcalcium), Cint, (Ptr{Int}, Ptr{Int}, Ref{qqbar}), p, q, a))
      throw(DomainError(a, "function value is not algebraic"))
   end
   return qqbar(p[1]) // q[1]
end



###############################################################################
#
#   Guessing
#
###############################################################################

function guess(R::CalciumQQBarField, x::acb, maxdeg::Int, maxbits::Int=0)
   prec = precision(parent(x))
   if maxbits <= 0
      maxbits = prec
   end
   res = qqbar()
   found = Bool(ccall((:qqbar_guess, libcalcium), Cint, (Ref{qqbar}, Ref{acb}, Int, Int, Int, Int), res, x, maxdeg, maxbits, 0, prec))
   if !found
      error("No suitable algebraic number found")
   end
   return res
end

function guess(R::CalciumQQBarField, x::arb, maxdeg::Int, maxbits::Int=0)
   CC = AcbField(precision(parent(x)))
   return guess(R, CC(x), maxdeg, maxbits)
end

###############################################################################
#
#   Conversions
#
###############################################################################

function (R::ArbField)(a::qqbar)
   prec = precision(R)
   z = R()
   ccall((:qqbar_get_arb, libcalcium), Nothing, (Ref{arb}, Ref{qqbar}, Int), z, a, prec)
   !isfinite(z) && throw(DomainError(a, "nonreal algebraic number"))
   return z
end

function (R::AcbField)(a::qqbar)
   prec = precision(R)
   z = R()
   ccall((:qqbar_get_acb, libcalcium), Nothing, (Ref{acb}, Ref{qqbar}, Int), z, a, prec)
   return z
end

# todo: provide qqbar_get_fmpq, qqbar_get_fmpz in C
function fmpq(a::qqbar)
   !isrational(a) && throw(DomainError(a), "nonrational algebraic number")
   p = fmpz()
   q = fmpz()
   ccall((:fmpz_poly_get_coeff_fmpz, libflint), Nothing, (Ref{fmpz}, Ref{qqbar}, Int), p, a, 0)
   ccall((:fmpz_poly_get_coeff_fmpz, libflint), Nothing, (Ref{fmpz}, Ref{qqbar}, Int), q, a, 1)
   ccall((:fmpz_neg, libflint), Nothing, (Ref{fmpz}, Ref{fmpz}), p, p)
   return p // q
end

function fmpz(a::qqbar)
   !isinteger(a) && throw(DomainError(a), "noninteger algebraic number")
   z = fmpz()
   ccall((:fmpz_poly_get_coeff_fmpz, libflint), Nothing, (Ref{fmpz}, Ref{qqbar}, Int), z, a, 0)
   ccall((:fmpz_neg, libflint), Nothing, (Ref{fmpz}, Ref{fmpz}), z, z)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(z::qqbar)
   ccall((:qqbar_zero, libcalcium), Nothing, (Ref{qqbar},), z)
   return z
end

function mul!(z::qqbar, x::qqbar, y::qqbar)
   ccall((:qqbar_mul, libcalcium), Nothing,
                (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, x, y)
   return z
end

function addeq!(z::qqbar, x::qqbar)
   ccall((:qqbar_add, libcalcium), Nothing,
                (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, z, x)
   return z
end

function add!(z::qqbar, x::qqbar, y::qqbar)
   ccall((:qqbar_add, libcalcium), Nothing,
                (Ref{qqbar}, Ref{qqbar}, Ref{qqbar}), z, x, y)
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

(a::CalciumQQBarField)(b::qqbar) = b

