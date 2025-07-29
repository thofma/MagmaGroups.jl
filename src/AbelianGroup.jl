mutable struct MagmaFinGenAbGroup <: MagmaGroup
	x::MagmaObject
  map

  MagmaFinGenAbGroup(x::MagmaObject) = new(x)

  MagmaFinGenAbGroup(x::MagmaObject, y) = new(x, y)
end

mutable struct MagmaFinGenAbGroupAutGrp <: MagmaGroup
	A::MagmaFinGenAbGroup
	x::MagmaObject
end

mutable struct MagmaMap{S, T}
	D::S
	C::T
	x::MagmaObject
end

domain(f::MagmaMap) = f.D

codomain(f::MagmaMap) = f.C

function is_bijective(f::MagmaMap)
  @assert domain(f) === codomain(f)
  return magconvert(Int, mag"#Kernel($(f.x))") == 1
end

function image(f::MagmaMap, x)
	return MagmaGroupElem(codomain(f), mag"$(data(x))@$(f.x)")
end

function preimage(f::MagmaMap, x)
	return MagmaGroupElem(domain(f), mag"$(data(x))@@$(f.x)")
end

data(G::MagmaFinGenAbGroup) = G.x

data(G::MagmaFinGenAbGroupAutGrp) = G.x

function abelian_group(::Type{Magma}, v::Vector{<:IntegerUnion})
	x = mag"AbelianGroup$(_prepare_for_magma.(v))"
	return MagmaFinGenAbGroup(x)
end

function abelian_group(::Type{Magma}, r::ZZMatrix)
  n = ncols(r)
  rr = BigInt.(Matrix(r))
  A = mag"FreeAbelianGroup($n)"
  rels = [ magcoerce(A, _prepare_for_magma(rr[i, :])) for i in 1:nrows(r)]
  f = mag"B, f := quo<$A | $rels>; return f"
  return MagmaFinGenAbGroup(mag"Codomain($f)", f)
end

function automorphism_group(A::MagmaFinGenAbGroup)
	AG = mag"AutomorphismGroup($(data(A)))"
	return MagmaFinGenAbGroupAutGrp(A, AG)
end

function (A::MagmaFinGenAbGroup)(x::Vector)
	xx = _prepare_for_magma(x)
  if !isdefined(A, :map)
    return MagmaGroupElem(A, mag"$(data(A))!$(xx)")
  else
    MagmaGroupElem(A, mag"(Domain($(A.map))!$(xx))@$(A.map)")
  end
end

function hom(A::MagmaGroup, B::MagmaGroup, x::Vector)
	@assert all(y -> parent(y) === B, x)
	xdata = data.(x)
	Ad = data(A)
	Bd = data(B)
	y = mag"y := $(Ad); return y";
	h = mag"hom<$y -> $(Bd) | $(xdata)>"
	return MagmaMap(A, B, h)
end 

function kernel(f::MagmaMap)
  _K = magf.Kernel(f.x)
  if domain(f) isa MagmaMatrixGroup
    K = typeof(domain(f))(domain(f).F, _K)
  else
    K = typeof(domain(f))(_K)
  end
  KtoD = MagmaMap(K, domain(f), mag"func< x | x>")
  return K, KtoD
end

function (A::MagmaFinGenAbGroupAutGrp)(h::MagmaMap)
	@assert domain(h) === A.A === codomain(h)
	x = mag"$(A.x)!$(h.x)"
	return MagmaGroupElem(A, x)
end

function Base.show(io::IO, a::MagmaGroupElem{MagmaFinGenAbGroup})
  A = parent(a)
  if !isdefined(A, :map)
    print(io, magconvert(Vector, mag"Eltseq($(data(a)))"))
  else
    print(io, magconvert(Vector, mag"Eltseq($(data(a))@@$(A.map))"))
  end
end

function _eltseq(a::MagmaGroupElem{MagmaFinGenAbGroup})
  A = parent(a)
  if !isdefined(A, :map)
    return ZZ.(magconvert(Vector, mag"Eltseq($(data(a)))"))
  else
    return ZZ.(magconvert(Vector, mag"Eltseq($(data(a))@@$(A.map))"))
  end
end

function isomorphism(::Type{PermGroup}, A::MagmaFinGenAbGroupAutGrp)
	f, P, _ = magf3.PermutationRepresentation(data(A))
	return MagmaMap(A, MagmaPermGroup(P), f)
end

function id(A::MagmaFinGenAbGroup)
  return MagmaGroupElem(A, mag"Identity($(data(A)))")
end

elem_type(::Type{MagmaFinGenAbGroup}) = MagmaGroupElem{MagmaFinGenAbGroup}

function Base.getproperty(a::MagmaGroupElem{MagmaFinGenAbGroup}, sym::Symbol)
  if sym === :coeff
    return matrix(ZZ, [_eltseq(a)])
  else
    return getfield(a, sym)
  end
end

function hom(a::MagmaGroupElem{MagmaFinGenAbGroupAutGrp})
  A = parent(a)
  G = A.A # underlying group
  return MagmaMap(G, G, mag"Parent(IdentityHomomorphism($(data(G))))!$(a.x)")
end

function MagmaGroup(A::FinGenAbGroup)
  S, StoA = snf(A)
  e = BigInt.(S.snf)
  MA = abelian_group(Magma, e)

  AtoMA = function(x::FinGenAbGroupElem)
    @assert parent(x) === A
    co = BigInt.(Hecke._eltseq(preimage(StoA, x).coeff))
    return MagmaGroupElem(MA, mag"$(data(MA))!$(co)")
  end

  MAtoA = function(y)
    a = StoA(S(magconvert(Vector{BigInt}, mag"Eltseq($(data(y)))")))
    @assert parent(a) === A
    return a
  end

  for i in 1:10
    x = rand(A, 10)
    y = rand(A, 10)
    @assert MAtoA(AtoMA(x)) == x
    @assert AtoMA(x + y) == AtoMA(x) + AtoMA(y)
  end

  return MA, AtoMA, MAtoA
end

function Oscar.isomorphism(::Type{FinGenAbGroup}, A::MagmaFinGenAbGroup)
  Arel = mag"RelationMatrix($(data(A)))"
  fl = magb.IsDiagonal(Arel)
  @assert fl
  dia = magconvert(Vector, mag"Diagonal($(Arel))")
  B = abelian_group(dia)

  BtoA = function(x::FinGenAbGroupElem)
    @assert parent(x) === B
    co = BigInt.(Hecke._eltseq(x.coeff))
    return MagmaGroupElem(A, mag"$(data(A))!$(co)")
  end

  AtoB = function(y)
    a = B(magconvert(Vector{BigInt}, mag"Eltseq($(data(y)))"))
    @assert parent(a) === B
    return a
  end

  for i in 1:10
    x = rand(B, 10)
    y = rand(B, 10)
    @assert AtoB(BtoA(x)) == x
    @assert BtoA(x + y) == BtoA(x) + BtoA(y)
  end

  return B, AtoB, BtoA
end

function Base.:(+)(x::MagmaGroupElem, y::MagmaGroupElem)
  @assert parent(x) === parent(y)
  return MagmaGroupElem(parent(x), mag"$(data(x)) + $(data(y))")
end

function Base.:(==)(x::MagmaGroupElem, y::MagmaGroupElem)
  return data(x) == data(y)
end

function is_zero(a::MagmaGroupElem)
  return magb.IsIdentity(data(a))
end
