###############################################################################
#
#   PariTypes.jl : Pari types
#
###############################################################################

###############################################################################
#
#   PariIntegerRing / pari_int
#
###############################################################################

type PariIntegerRing <: PariRing
end

type pari_int <: RingElem
   d::Ptr{Int}

   function pari_int(s::Int)
      g = new(ccall((:pari_malloc, :libpari), Ptr{Int}, 
                    (Int,), s*BITS_IN_WORD))
      finalizer(g, _pari_int_clear_fn)
      return g
   end
end

function _pari_int_clear_fn(g::pari_int)
   ccall((:pari_free, :libpari), Void, (Ptr{Uint},), g.d)
end

###############################################################################
#
#   PariRationalField / pari_rat
#
###############################################################################

type PariRationalField <: PariField
end

type pari_rat <: RingElem
   d::Ptr{Int}

   function pari_rat(s::Int)
      g = new(ccall((:pari_malloc, :libpari), Ptr{Int}, 
                    (Int,), s*BITS_IN_WORD))
      finalizer(g, _pari_rat_clear_fn)
      return g
   end
end

_pari_rat_clear_fn(g::pari_rat) = ccall((:pari_free, :libpari), Void, 
                                        (Ptr{Uint},), g.d)

###############################################################################
#
#   PariVector / pari_vec
#
###############################################################################

type PariVector{T <: Union(PariRing, PariSet)}
   base_ring::Union(PariRing, PariSet)
end

type pari_vec{T <: Union(PariRing, PariSet)}
   data::Ptr{Int}
   parent::PariVector{T}

   function pari_vec{R <: Union(PariRing, PariSet)}(v::Ptr{Int}, par::R)
      r = new(gclone(v), PariVector{T}(par))
      finalizer(r, _pari_vec_unclone)
      return r
   end
end

_pari_vec_unclone(a::pari_vec) = gunclone(a.data)

###############################################################################
#
#   PariPolyRing / pari_poly
#
###############################################################################

PariPolyID = ObjectIdDict()

type PariPolyRing{T <: PariRing} <: PariRing
   base_ring::PariRing
   pol_0::Ptr{Int}
   S::Symbol

   function PariPolyRing(R::PariRing, s::Symbol)
      z = ccall((:pari_malloc, :libpari), Ptr{Int}, (Int,), 2*sizeof(Int))
      unsafe_store!(z, evaltyp(t_POL) | 2, 1) 
      unsafe_store!(z, evalsigne(0) | evalvarn(0), 2)
      try
         return PariPolyID[R, s]
      catch
         r = PariPolyID[R, s] = new(R, z, s)
         finalizer(r, _pari_poly_zero_clear_fn)
         return r
      end
      
   end
end

function _pari_poly_zero_clear_fn(p::PariPolyRing)
   ccall((:pari_free, :libpari), Void, (Ptr{Int},), p.pol_0)
end

type pari_poly{T <: PariRing} <: PolyElem{T}
   d::Ptr{Int}
   parent::PariPolyRing{T}

   function pari_poly(data::Ptr{Int})
      g = new(gclone(data))
      finalizer(g, _pari_poly_unclone)
      return g
   end

   function pari_poly(s::Int)
      g = new(ccall((:pari_malloc, :libpari), Ptr{Int}, (Int,), s*sizeof(Int)))
      finalizer(g, _pari_poly_clear_fn)
      return g
   end
end

function _pari_poly_clear_fn(g::pari_poly)
   ccall((:pari_free, :libpari), Void, (Ptr{Uint},), g.d)
end

_pari_poly_unclone(g::pari_poly) = gunclone(g.d)

###############################################################################
#
#   PariPolModRing / pari_polmod
#
###############################################################################

PariPolModID = Dict{Tuple{DataType, Symbol}, PariRing}()

type PariPolModRing{S <: PariRing} <: PariRing
   T::Symbol

   function PariPolModRing(t::Symbol)
      try
         return PariPolModID[S, t]
      catch
         return PariPolModID[S, t] = new(t)
      end
   end
end

type pari_polmod{S <: PariRing} <: PolyElem{S}
   data::Ptr{Int}
   parent::PariPolModRing{S}

   function pari_polmod(data::Ptr{Int})
      r = new(gclone(data))
      finalizer(r, _pari_polmod_unclone)
      return r
   end
end

_pari_polmod_unclone(a::pari_polmod) = gunclone(a.data)

###############################################################################
#
#   PariNumberField
#
###############################################################################

PariNumberFieldID = Dict{fmpq_poly, PariRing}()

type PariNumberField <: PariRing
   data::Ptr{Int}
   nf::NfNumberField
   
   function PariNumberField(nf::NfNumberField)
      try
         return PariNumberFieldID[nf.pol]
      catch
         av = unsafe_load(avma, 1)
         p = pari(nf.pol)
         d = gclone(ccall((:nfinit, :libpari), Ptr{Int}, 
                           (Ptr{Int}, Int), p.d, 5))
         unsafe_store!(avma, av, 1)
         ord = new(d, nf)
         finalizer(ord, _pari_nf_unclone)
         return PariNumberFieldID[nf.pol] = ord
      end
   end
end

_pari_nf_unclone(a::PariNumberField) = gunclone(a.data)

###############################################################################
#
#   PariMaximalOrder / PariMaximalOrderElem
#
###############################################################################

type PariMaximalOrder <: PariRing
   pari_nf::PariNumberField
end

type PariMaximalOrderElem <: RingElem
   data::Ptr{Int}
   parent::PariMaximalOrder

   function PariMaximalOrderElem(a::Ptr{Int}, par::PariMaximalOrder)
      r = new(gclone(a), par)
      finalizer(r, _pari_maximal_order_elem_clear_fn)
      return r
   end
end

_pari_maximal_order_elem_clear_fn(a::PariMaximalOrderElem) = gunclone(a.data)

###############################################################################
#
#   PariIdealSet / PariIdeal
#
###############################################################################

PariIdealSetID = ObjectIdDict()

type PariIdealSet <: PariSet
   order::PariMaximalOrder

   function PariIdealSet(ord::PariMaximalOrder)
      return try
         PariIdealSetID[ord]
      catch
         PariIdealSetID[ord] = new(ord)
      end
   end
end

type PariIdeal <: PariSet
   ideal::Ptr{Int}
   parent::PariIdealSet

   function PariIdeal(a::Ptr{Int}, par::PariIdealSet)
      r = new(gclone(a), par)
      finalizer(r, _pari_ideal_clear_fn)
      return r
   end
end

_pari_ideal_clear_fn(a::PariIdeal) = gunclone(a.ideal)