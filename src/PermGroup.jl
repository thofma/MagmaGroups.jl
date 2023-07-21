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
