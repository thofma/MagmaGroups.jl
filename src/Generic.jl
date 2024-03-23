import Oscar: order

order(x::MagmaGroup) = ZZRingElem(magconvert(BigInt, mag"Order($(data(x)))"))

ngens(x::MagmaGroup) = magconvert(Int, mag"Ngens($(data(x)))")

mutable struct MagmaGroupElem{S}
	G::S
	x::MagmaObject
end

data(g::MagmaGroupElem) = g.x

parent(x::MagmaGroupElem) = x.G

function gen(A::MagmaGroup, i::Int)
	return MagmaGroupElem(A, mag"$(data(A)).$(i)")
end

function sub(G::MagmaGroup, x::Vector{<: MagmaGroupElem})
	xx = data.(x)
	_f = mag"H, HtoG := sub< $(data(G)) | $(xx)>; return HtoG"
  if G isa MagmaMatrixGroup
    H = typeof(G)(G.F, mag"Domain($(_f))")
  else
    H = typeof(G)(mag"Domain($(_f))")
  end
  HtoG = MagmaMap(H, G, _f)
  return H, HtoG
end

function rand(G::MagmaGroup)
  return MagmaGroupElem(G, mag"Random($(data(G)))")
end

function maximal_abelian_quotient(G::MagmaGroup)
  _Q, _f = magf2.AbelianQuotient(data(G))
  Q = MagmaGrpAbFinGen(_Q)
  return Q, MagmaMap(G, Q, _f)
end

function Base.collect(x::T) where {T <: MagmaGroup}
  _r =  mag"[y : y in $(x.x)]"
  return MagmaGroupElem{T}[MagmaGroupElem(x, mag"$(_r)[$i]") for i in 1:Int(order(x))]
end

function double_coset_representatives(G::MagmaGroup, H::MagmaGroup, K::MagmaGroup)
  x = mag"DoubleCosetRepresentatives($(data(G)), $(data(H)), $(data(K)))"
  yy = MagmaGroupElem{typeof(G)}[]
  k = magconvert(Int, mag"#$(x)");
  for i in 1:k
    push!(yy, MagmaGroupElem(G, mag"$(x)[$i]"))
  end
  return yy
end
