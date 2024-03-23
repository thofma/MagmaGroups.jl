mutable struct MagmaMatrixGroup{S, T} <: MagmaGroup
  F
  x::MagmaObject
end

data(G::MagmaMatrixGroup) = G.x

function Oscar.Hecke._eltseq(x::MagmaGroupElem{<:MagmaMatrixGroup})
  v = mag"[Integers()!x : x in Eltseq($(data(x)))]" 
  n = magi.Nrows(data(x))
  F = parent(x).F
  return F.(magconvert(Vector, v))
end

function (G::MagmaMatrixGroup)(y::MatElem)
  @assert base_ring(y) === G.F
  _m = oscar_matrix_to_magma_matrix(mag"BaseRing($(data(G)))", y)
  return MagmaGroupElem(G, mag"$(data(G))!$(_m)")
end

function oscar_gf_to_magma_gf(F)
  p = Int(characteristic(F))
  return magf.GF(p)
end

function matrix_group(::Type{Magma}, V::Vector{<:MatElem})
  F = base_ring(V[1])
  _F = oscar_gf_to_magma_gf(F)
  @info "Translating matrices to Magma"
  vv = magseq()
  for m in V
    magp1.Append(vv, oscar_matrix_to_magma_matrix(_F, m))
  end
  n = nrows(V[1])
  @info "Creating group"
  S = mag"MatrixGroup<$n, $(_F) | $vv>";
  return MagmaMatrixGroup{elem_type(F), dense_matrix_type(F)}(F, S)
end

function oscar_matrix_to_magma_matrix(F, m)
  n = nrows(m)
  _m = map(x -> BigInt(lift(x)), Oscar.Hecke._eltseq(m))
  v = magf.Matrix(F, n, n, _m)
  #v = mag"Matrix($F, $n, $n, $(_m))" 
  return v
end

function magma_matrix_group(v::Vector)
end

function magma_matrix_to_oscar_matrix(F, m)
  v = mag"[Integers()!x : x in Eltseq($m)]" 
  n = magi.Nrows(m)
  return matrix(F, n, n, ZZ.(magconvert(Vector, v)))
end


function magma_matrix_to_matrix_group(S, m)
  return mag"$(S)!$(m)"
end


