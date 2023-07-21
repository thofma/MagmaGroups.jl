import Oscar: order

order(x::MagmaGroup) = magconvert(Int, mag"Order($(data(x)))")

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
	mag"H, HtoG := sub< $(data(G)) | $(xx)>; return HtoG"
end

