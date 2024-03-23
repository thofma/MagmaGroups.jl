mutable struct MagmaPermGroup <: MagmaGroup
	x::MagmaObject
end

data(G::MagmaPermGroup) = G.x

function permutation_group(::Type{Magma}, n::Int, x::Vector{PermGroupElem})
	lists = [[x[i](j) for j in 1:n] for i in 1:length(x)]
	return MagmaPermGroup(mag"PermutationGroup< $n | $lists>")
end

function (G::MagmaPermGroup)(x::Vector{Int})
	return MagmaGroupElem(G, mag"$(data(G))!$x")
end

degree(G::MagmaPermGroup) = magconvert(Int, mag"Degree($(data(G)))")

function inner_direct_product(Gs::Vector{MagmaPermGroup}; morphisms = true)
  @assert morphisms
  Gss = data.(Gs)
  _P, _injs, _projs = magf3.DirectProduct(Gss)
  P = MagmaPermGroup(_P)
  injs = [ MagmaMap(Gs[i], P, _injs[i]) for i in 1:length(Gs) ]
  projs = [ MagmaMap(P, Gs[i], _projs[i]) for i in 1:length(Gs) ]
  return P, injs, projs
end

function one(G::MagmaPermGroup)
  return MagmaGroupElem(G, mag"Identity($(data(G)))")
end

function Base.:(*)(x::MagmaGroupElem, y::MagmaGroupElem)
  @assert parent(x) === parent(y)
  return MagmaGroupElem(parent(x), mag"$(data(x)) * $(data(y))")
end
