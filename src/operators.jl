function hamiltonian(L::Lattice, lH::AbstractOperator{B,B}, tH::Vector{Tuple{Vararg{O,2}}}) where {B<:Basis,O<:AbstractOperator{B,B}}
    @assert L.lbasis == lH.basis_l "Local Hamiltonian's basis incompatible with lattice's local basis."
    H = begin
            if (typeof(lH) <: SparseOperator) && (typeof(lH) <: SparseOperator)
                SparseOperator(L.gbasis);
            else
                DenseOperator(L.gbasis);
            end
        end
    for e in edges(L)
        for i in 1:length(tH)
            H.data .+= embed(L.gbasis,[e.src, e.dst],collect(tH[i])).data;
        end
    end
    H.data .= H.data + H.data';
    @inbounds for v in vertices(L)
        H.data .+= embed(L.gbasis, [v], [lH]).data;
    end
    return H;
end

function dissipators(L::Lattice, J::Vector{O}) where {B<:Basis,O<:AbstractOperator{B,B}}
    @assert L.lbasis == first(J).basis_l "Local jump operators' basis incompatible with lattice's local basis."
    return union([[embed(L.gbasis, [v], [J[i]]) for i in 1:length(J)] for v in vertices(L)]...)
end
