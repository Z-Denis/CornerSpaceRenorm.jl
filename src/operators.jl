"""
    hamiltonian(L, lH, tH)

Construct a Hamiltonian from a `Lattice` and operators.
# Arguments
* `L`: Lattice.
* `lH`: local Hamiltonian.
* `tH`: `Vector` of `Tuple` of two local operators defining the tunneling.
"""
function hamiltonian(L::Lattice, lH::AbstractOperator{B,B}, tH::Vector{Tuple{Vararg{O,2}}}) where {B<:Basis,O<:AbstractOperator{B,B}}
    gbasis = CompositeBasis([lH.basis_l for i in 1:nv(L)]...);
    H = begin
            if (typeof(lH) <: SparseOperator) && (typeof(lH) <: SparseOperator)
                SparseOperator(gbasis);
            else
                DenseOperator(gbasis);
            end
        end
    for e in edges(L)
        for i in 1:length(tH)
            H.data .+= embed(gbasis,[e.src, e.dst],collect(tH[i])).data;
        end
    end
    H.data .= H.data + H.data';
    @inbounds for v in vertices(L)
        H.data .+= embed(gbasis, [v], [lH]).data;
    end
    return H;
end

#hamiltonian(s::System) = hamiltonian(s.lattice, s.lH, [s.tH])

"""
    dissipators(L, J)

Construct the jump operators of a system from a `Lattice` and local jump operators.
# Arguments
* `L`: Lattice.
* `J`: Array of on-site local jump operators.
"""
function dissipators(L::Lattice, J::Vector{O}) where {B<:Basis,O<:AbstractOperator{B,B}}
    gbasis = CompositeBasis([first(J).basis_l for i in 1:nv(L)]...);
    return Base.union([[embed(gbasis, [v], [J[i]]) for i in 1:length(J)] for v in vertices(L)]...)
end

#dissipators(s::System) = dissipators(s.lattice, s.J)
