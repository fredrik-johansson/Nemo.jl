###############################################################################
#
#   arb_abs_series.jl : Power series over arb reals (using arb_poly)
#
###############################################################################

export arb_abs_series, ArbAbsSeriesRing, PowerSeriesRing,
       exp, log, tan, tanh, sin, sinh, cos, cosh, asin, asinh, atan,
       atanh, acos, acosh, sqrt, rsqrt

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

function O(a::arb_abs_series)
   prec = length(a) - 1
   prec < 0 && throw(DomainError())
   z = parent(a)()
   z.prec = prec
   z.parent = parent(a)
   return z
end

elem_type(::ArbAbsSeriesRing) = arb_abs_series

parent_type(::Type{arb_abs_series}) = ArbAbsSeriesRing

base_ring(R::ArbAbsSeriesRing) = R.base_ring

var(a::ArbAbsSeriesRing) = a.S

###############################################################################
#
#   Basic manipulation
#
###############################################################################    
   
max_precision(R::ArbAbsSeriesRing) = R.prec_max

function normalise(a::arb_abs_series, len::Int)
   if len > 0
      c = base_ring(parent(a))()
      ccall((:arb_poly_get_coeff_arb, :libarb), Void, 
         (Ptr{arb}, Ptr{arb_abs_series}, Int), &c, &a, len - 1)
   end
   while len > 0 && iszero(c)
      len -= 1
      if len > 0
         ccall((:arb_poly_get_coeff_arb, :libarb), Void, 
            (Ptr{arb}, Ptr{arb_abs_series}, Int), &c, &a, len - 1)
      end
   end

   return len
end

function coeff(x::arb_abs_series, n::Int)
   if n < 0
      return base_ring(parent(x))()
   end
   z = base_ring(parent(x))()
   ccall((:arb_poly_get_coeff_arb, :libarb), Void, 
         (Ptr{arb}, Ptr{arb_abs_series}, Int), &z, &x, n)
   return z
end

function length(x::arb_abs_series)
   return ccall((:arb_poly_length, :libarb), Int, (Ptr{arb_abs_series},), &x)
end

precision(x::arb_abs_series) = x.prec

zero(R::ArbAbsSeriesRing) = R(0)

one(R::ArbAbsSeriesRing) = R(1)

function gen(R::ArbAbsSeriesRing)
   r = base_ring(R)
   z = arb_abs_series([r(0), r(1)], 2, max_precision(R))
   z.parent = R
   return z
end

function deepcopy(a::arb_abs_series)
   z = arb_abs_series(a)
   z.prec = a.prec
   z.parent = parent(a)
   return z
end

function isgen(a::arb_abs_series)
   return precision(a) == 0 || (precision(a) == 1 && iszero(a)) ||
      ccall((:arb_poly_is_x, :libarb), Bool, (Ptr{arb_abs_series}, ), &a)
end

iszero(a::arb_abs_series) = length(a) == 0

isunit(a::arb_abs_series) = valuation(a) == 0 && isunit(coeff(a, 0))

function isone(a::arb_abs_series)
   return precision(a) == 0 ||
      ccall((:arb_poly_is_one, :libarb), Bool, (Ptr{arb_abs_series}, ), &a)
end

function valuation(a::arb_abs_series)
   v = ccall((:arb_poly_valuation, :libarb), Int, (Ptr{arb_abs_series}, ), &a)
   if v == -1
      return precision(a)
   else
      return v
   end
end

###############################################################################
#
#   AbstractString{} I/O
#
###############################################################################

function show(io::IO, x::arb_abs_series)
   len = length(x)
   if len == 0
      print(io, zero(base_ring(x)))
   else
      coeff_printed = false
      for i = 0:len - 1
         c = coeff(x, i)
         if !iszero(c)
            if coeff_printed
               print(io, "+")
            end
            if i != 0
               if !isone(c)
                  print(io, "(")
                  print(io, c)
                  print(io, ")")
                  if i != 0
                     print(io, "*")
                  end
               end
               print(io, string(var(parent(x))))
               if i != 1
                  print(io, "^")
                  print(io, i)
               end
            else
               print(io, c)
            end
            coeff_printed = true
         end
      end
   end
   print(io, "+O(", string(var(parent(x))), "^", precision(x), ")")
end


function show(io::IO, a::ArbAbsSeriesRing)
   print(io, "Univariate power series ring in ", var(a), " over ")
   show(io, base_ring(a))
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(x::arb_abs_series)
   z = parent(x)()
   ccall((:arb_poly_neg, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}), 
               &z, &x)
   z.prec = x.prec
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::arb_abs_series, b::arb_abs_series)
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
         
   p = min(a.prec, b.prec)
   z = parent(a)()
   z.prec = p

   bitprec = prec(parent(a))

   ccall((:arb_poly_add_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, &b, p, bitprec)
   return z
end

function -(a::arb_abs_series, b::arb_abs_series)
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
         
   p = min(a.prec, b.prec)
   z = parent(a)()
   z.prec = p

   bitprec = prec(parent(a))

   ccall((:arb_poly_sub_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, &b, p, bitprec)
   return z
end

function *(a::arb_abs_series, b::arb_abs_series)

   check_parent(a, b)
   lena = length(a)
   lenb = length(b)

   aval = valuation(a)
   bval = valuation(b)

   p = min(a.prec + bval, b.prec + aval)
   p = min(p, max_precision(parent(a)))
   
   lena = min(lena, p)
   lenb = min(lenb, p)
   
   z = parent(a)()
   z.prec = p
      
   if lena == 0 || lenb == 0
      return z
   end

   lenz = min(lena + lenb - 1, p)

   bitprec = prec(parent(a))

   ccall((:arb_poly_mullow, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, &b, lenz, bitprec)

   return z
end


###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function divexact(x::arb_abs_series, y::Union{Int,fmpz,fmpq,arb})
   y == 0 && throw(DivideError())
   y = base_ring(x)(y)
   z = parent(x)()
   z.prec = x.prec
   bitprec = prec(parent(x))
   ccall((:arb_poly_scalar_div, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb}, Int), 
               &z, &x, y, bitprec)
   return z
end

*(x::arb_abs_series, y::Union{Int,fmpz,fmpq,arb}) = y*x

###############################################################################
#
#   Shifting
#
###############################################################################

function shift_left(x::arb_abs_series, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   z = parent(x)()
   z.prec = x.prec + len
   ccall((:arb_poly_shift_left, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int), 
               &z, &x, len)
   return z
end

function shift_right(x::arb_abs_series, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   z = parent(x)()
   if len >= xlen
      z.prec = max(0, x.prec - len)
   else
      z.prec = x.prec - len
      ccall((:arb_poly_shift_right, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int), 
               &z, &x, len)
   end
   return z
end

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate(x::arb_abs_series, prec::Int)
   prec < 0 && throw(DomainError())
   if x.prec <= prec
      return x
   end
   z = parent(x)()
   z.prec = prec
   ccall((:arb_poly_set_trunc, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int), 
               &z, &x, prec)
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::arb_abs_series, b::arb)
   p = a.prec
   z = parent(a)()
   z.prec = p
   bitprec = prec(parent(a))
   ccall((:arb_poly_pow_arb_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb}, Int, Int), 
               &z, &a, &b, p, bitprec)
   return z
end

function ^(a::arb_abs_series, b::arb_abs_series)
   check_parent(a, b)
   p = min(a.prec, b.prec)
   z = parent(a)()
   z.prec = p
   bitprec = prec(parent(a))
   ccall((:arb_poly_pow_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, &b, p, bitprec)
   return z
end

function ^(a::arb_abs_series, b::Int)
   b < 0 && return a ^ base_ring(parent(a))(b)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_pow_ui_trunc_binexp, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, UInt, Int, Int), 
               &z, &a, b, z.prec, bitprec)
   set_prec!(z, precision(a))
   return z
end

# todo: more coercions

###############################################################################
#
#   Comparison
#
###############################################################################

# todo ==, !=, isequal

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::arb_abs_series, y::arb_abs_series)
   check_parent(x, y)
   y == 0 && throw(DivideError())
   v2 = valuation(y)
   v1 = valuation(x)
   if v2 != 0
      if v1 >= v2
         x = shift_right(x, v2)
         y = shift_right(y, v2)
      end
   end
   !isunit(y) && error("Unable to invert power series")
   prec = min(x.prec, y.prec - v2 + v1)
   z = parent(x)()
   z.prec = prec
   bitprec = prec(parent(x))
   ccall((:arb_poly_div_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &x, &y, prec, bitprec)
   return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::arb_abs_series, y::Union{Int,Integer,fmpz,fmpq,arb})
   y == 0 && throw(DivideError())
   y = base_ring(x)(y)
   z = parent(x)()
   z.prec = x.prec
   bitprec = prec(parent(x))
   ccall((:arb_poly_scalar_div, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb}, Int), 
               &z, &x, y, bitprec)
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::arb_abs_series)
   a == 0 && throw(DivideError())
   ainv = parent(a)()
   ainv.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_inv_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &ainv, &a, a.prec, bitprec)
   return ainv
end

###############################################################################
#
#   Special functions
#
###############################################################################

function exp(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_exp_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function log(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_log_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function tan(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_tan_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function tanh(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_tanh_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function sin(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_sin_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function sinh(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_sinh_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function cos(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_cos_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function cosh(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_cosh_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function asin(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_asin_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function asinh(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_asinh_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function atan(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_atan_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function atanh(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_atanh_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function sqrt(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_sqrt_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

function rsqrt(a::arb_abs_series)
   z = parent(a)()
   z.prec = a.prec
   bitprec = prec(parent(a))
   ccall((:arb_poly_rsqrt_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, a.prec, bitprec)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function setcoeff!(z::arb_abs_series, n::Int, x::arb)
   ccall((:arb_poly_set_coeff_arb, :libarb), Void, 
                (Ptr{arb_abs_series}, Int, Ptr{arb}), 
               &z, n, &x)
end

function mul!(z::arb_abs_series, a::arb_abs_series, b::arb_abs_series)
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   p = min(a.prec + bval, b.prec + aval)
   p = min(p, max_precision(parent(a)))

   lena = min(lena, p)
   lenb = min(lenb, p)
   
   lenz = min(lena + lenb - 1, p)
   if lenz < 0
      lenz = 0
   end

   bitprec = prec(parent(x))

   z.prec = p
   ccall((:arb_poly_mullow, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &z, &a, &b, lenz, bitprec)
end

function addeq!(a::arb_abs_series, b::arb_abs_series)
   lena = length(a)
   lenb = length(b)
         
   p = min(a.prec, b.prec)
 
   lena = min(lena, p)
   lenb = min(lenb, p)

   lenz = max(lena, lenb)

   bitprec = prec(parent(x))

   a.prec = p
   ccall((:arb_poly_add_series, :libarb), Void, 
                (Ptr{arb_abs_series}, Ptr{arb_abs_series}, Ptr{arb_abs_series}, Int, Int), 
               &a, &a, &b, lenz, bitprec)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

Base.promote_rule{T <: Integer}(::Type{arb_abs_series}, ::Type{T}) = arb_abs_series

Base.promote_rule(::Type{arb_abs_series}, ::Type{fmpz}) = arb_abs_series

Base.promote_rule(::Type{arb_abs_series}, ::Type{fmpq}) = arb_abs_series

Base.promote_rule(::Type{arb_abs_series}, ::Type{arb}) = arb_abs_series

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call(a::ArbAbsSeriesRing)
   z = arb_abs_series()
   z.prec = a.prec_max
   z.parent = a
   return z
end

function Base.call(a::ArbAbsSeriesRing, b::Union{Int,Integer,fmpz,fmpq,arb})
   if b == 0
      z = arb_abs_series()
      z.prec = a.prec_max
   else
      z = arb_abs_series([base_ring(a)(b)], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function Base.call(a::ArbAbsSeriesRing, b::arb_abs_series)
   parent(b) != a && error("Unable to coerce power series")
   return b
end

function Base.call(a::ArbAbsSeriesRing, b::Array{arb, 1}, len::Int, prec::Int)
   z = arb_abs_series(b, len, prec)
   z.parent = a
   return z
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

function PowerSeriesRing(R::ArbField, prec::Int, s::AbstractString{}; cached=true, model=:capped_relative)
   S = Symbol(s)

   if model == :capped_relative
      parent_obj = GenAbsSeriesRing(R, prec, S, cached)
   elseif model == :capped_absolute
      parent_obj = ArbAbsSeriesRing(R, prec, S)
   end

   return parent_obj, gen(parent_obj)
end

