module MagmaGroups

# Write your package code here.

using Oscar
using MagmaCall

import Oscar: IntegerUnion, automorphism_group, abelian_group, isomorphism, ngens, permutation_group, hom, domain, codomain, gen, image, preimage, parent, degree, inner_direct_product, id, elem_type, is_bijective, one, *, hom, sub, rand, matrix_group, maximal_abelian_quotient, +, kernel, is_zero

_prepare_for_magma(x::ZZRingElem) = BigInt(x)
_prepare_for_magma(x::Vector) = _prepare_for_magma.(x)
_prepare_for_magma(x) = x

abstract type Magma end

abstract type MagmaGroup end

include("Generic.jl")
include("PermGroup.jl")
include("AbelianGroup.jl")
include("MatrixGroup.jl")


export Magma

end
