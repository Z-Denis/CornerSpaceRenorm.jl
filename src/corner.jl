using DataStructures, LinearAlgebra

function maxk(a::AbstractVector, k::Integer)
    ix = partialsortperm(a, 1:k, rev=true)
    @inbounds @views collect(ix), collect(a[ix[1:k]])
end

function max_prod_pairs(a::AbstractVector, b::AbstractVector, k::Integer)
    ix_a, _a = maxk(a, min(k,length(a)))
    ix_b, _b = maxk(b, min(k,length(b)))
    heap_idcs = [CartesianIndex(i,1) for i in 1:length(_a)]
    h = MutableBinaryMaxHeap([_a[idx[1]] * _b[idx[2]] for idx in heap_idcs])
    handles = Array{Tuple{Int64,Int64},1}(undef, k)
    prod_pairs = zeros(eltype(a), k)
    i = 1
    while i <= k || length(h) == 0
        prod_pairs[i], heap_idx = top_with_handle(h);
        handles[i] = (ix_a[heap_idcs[heap_idx][1]], ix_b[heap_idcs[heap_idx][2]])
        pop!(h);
        if heap_idcs[heap_idx][2] < length(_b)
            push!(heap_idcs, CartesianIndex(heap_idcs[heap_idx][1], heap_idcs[heap_idx][2]+1))
            push!(h, _a[heap_idcs[heap_idx][1]] * _b[heap_idcs[heap_idx][2]+1])
        end
        i += 1
    end
    return handles, prod_pairs
end

function corner_subspace(ρA::DenseOperator{B1,B1}, ρB::DenseOperator{B2,B2}, k::Integer) where {B1<:Basis,B2<:Basis}
    bA = ρA.basis_l
    bB = ρB.basis_l
    ps_A, αs_A = eigen(ρA.data); # αs_A[:,i] : i-th eigenvector
    ps_B, αs_B = eigen(ρB.data);
    @inbounds ϕs_A = [Ket(bA, αs_A[:,i]) for i in 1:length(ps_A)]
    @inbounds ϕs_B = [Ket(bB, αs_B[:,i]) for i in 1:length(ps_B)]
    handles, prod_pairs = max_prod_pairs(real.(ps_A), real.(ps_B), k)
    bC = typeof(bA) <: CompositeBasis ? CompositeBasis(bA.bases...,bB.bases...) : CompositeBasis(bA,bB)
    return SubspaceBasis(bC, [ϕs_A[idcs[1]] ⊗ ϕs_B[idcs[2]] for idcs in handles]), prod_pairs
end
