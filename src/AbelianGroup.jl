mutable struct MagmaGrpAbFinGen <: MagmaGroup
	x::MagmaObject
end

mutable struct MagmaGrpAbFinGenAutGrp <: MagmaGroup
	A::MagmaGrpAbFinGen
	x::MagmaObject
end

mutable struct MagmaMap{S, T}
	D::S
	C::T
	x::MagmaObject
end

domain(f::MagmaMap) = f.D

codomain(f::MagmaMap) = f.C

function image(f::MagmaMap, x)
	return MagmaGroupElem(codomain(f), mag"$(data(x))@$(f.x)")
end

function preimage(f::MagmaMap, x)
	return MagmaGroupElem(domain(f), mag"$(data(x))@@$(f.x)")
end

data(G::MagmaGrpAbFinGen) = G.x

data(G::MagmaGrpAbFinGenAutGrp) = G.x

function abelian_group(::Type{Magma}, v::Vector{<:IntegerUnion})
	x = mag"AbelianGroup$(_prepare_for_magma.(v))"
	return MagmaGrpAbFinGen(x)
end

function automorphism_group(A::MagmaGrpAbFinGen)
	AG = mag"AutomorphismGroup($(data(A)))"
	return MagmaGrpAbFinGenAutGrp(A, AG)
end

function (A::MagmaGrpAbFinGen)(x::Vector)
	xx = _prepare_for_magma(x)
	return MagmaGroupElem(A, mag"$(data(A))!$(xx)")
end

function hom(A::MagmaGrpAbFinGen, B::MagmaGrpAbFinGen, x::Vector)
	@assert all(y -> parent(y) === B, x)
	xdata = data.(x)
	h = mag"hom<$(data(A)) -> $(data(B)) | $(xdata)>"
	return MagmaMap(A, B, h)
end 

function (A::MagmaGrpAbFinGenAutGrp)(h::MagmaMap)
	@assert domain(h) === A.A === codomain(h)
	x = mag"$(A.x)!$(h.x)"
	return MagmaGroupElem(A, x)
end

function isomorphism(::Type{PermGroup}, A::MagmaGrpAbFinGenAutGrp)
	f, P, _ = magf3.PermutationRepresentation(data(A))
	return MagmaMap(A, MagmaPermGroup(P), f)
end
