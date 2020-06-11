using DataStructures, LinearAlgebra

"""
    CornerSpaceRenorm.cgs(X)

Apply the classical Gram-Schmidt orthogonalisation procedure to a set of vectors.
Vectors are passed as columns of a dense matrix `X`.
"""
function cgs(X::Matrix{T}) where T <: Number
    m, n = size(X);
    V = copy(X);
    Q = zeros(T, size(X));

    @inbounds @views for j = 1:n
        for k = 1:j-1
            axpy!(-BLAS.dotc(m, Q[:,k], 1, V[:,j], 1), Q[:,k], V[:,j])
        end
        Q[:,j] .= V[:,j] ./ sqrt(BLAS.dotc(n, V[:,j], 1, V[:,j], 1));
    end

    return Q
end

"""
    CornerSpaceRenorm.cgs!(X)

In-place classical Gram-Schmidt orthogonalisation procedure.
See also: [`CornerSpaceRenorm.cgs`](@ref)
"""
function cgs!(X::Matrix{T}) where T <: Number
    m, n = size(X);
    V = copy(X);

    @inbounds @views for j = 1:n
        for k = 1:j-1
            axpy!(-BLAS.dotc(m, X[:,k], 1, V[:,j], 1), X[:,k], V[:,j])
        end
        X[:,j] .= V[:,j] ./ sqrt(BLAS.dotc(n, V[:,j], 1, V[:,j], 1));
    end

    return X
end

"""
    CornerSpaceRenorm.mgs(X)

Apply the modified Gram-Schmidt orthogonalisation procedure to a set of vectors.
Vectors are passed as columns of a dense matrix `X`.
"""
function mgs(X::Matrix{T}) where T <: Number
    m, n = size(X);
    V = copy(X);
    Q = copy(X);

    for j = 1:n
        @inbounds @views Q[:,j] .= V[:,j]./norm(V[:,j]);
        @inbounds @views for k = j+1:n
            axpy!(-BLAS.dotc(m, Q[:,j], 1, V[:,k], 1), Q[:,j], V[:,k])
        end
    end

    return Q
end

"""
    CornerSpaceRenorm.mgs!(X)

In-place modified Gram-Schmidt orthogonalisation procedure.
See also: [`CornerSpaceRenorm.mgs`](@ref)
"""
function mgs!(X::Matrix{T}) where T <: Number
    m, n = size(X);
    V = copy(X);

    for j = 1:n
        @inbounds @views X[:,j] .= V[:,j]./norm(V[:,j]);
        @inbounds @views for k = j+1:n
            axpy!(-BLAS.dotc(m, X[:,j], 1, V[:,k], 1), X[:,j], V[:,k])
        end
    end

    return X
end

"""
    CornerSpaceRenorm.maxk(a, k)

Find `k` largest elements of an array and return their indices.
# Arguments
* `a`: some array.
* `k`: number of largest elements to be found.
"""
function maxk(a::AbstractVector{T}, k::Int) where T<:Real
    ix = partialsortperm(a, 1:k, rev=true)
    @inbounds @views collect(ix), collect(a[ix[1:k]])
end

"""
    CornerSpaceRenorm.max_prod_pairs(a, b, k)

Find `k` pairs `(a_i,b_i)` with largest products from two arrays `a` and `b`.
# Arguments
* `a`: some array.
* `b`: some array.
* `k`: number of largest products to be found.
"""
function max_prod_pairs(a::AbstractVector{T1}, b::AbstractVector{T2}, k::Int) where {T1<:Real, T2<:Real}
    ix_a, _a = maxk(a, min(k,length(a)))
    ix_b, _b = maxk(b, min(k,length(b)))
    heap_idcs = [CartesianIndex(i,1) for i in 1:length(_a)]
    h = MutableBinaryMaxHeap([_a[idx[1]] * _b[idx[2]] for idx in heap_idcs])
    handles = Array{Tuple{Int64,Int64},1}(undef, k)
    prod_pairs = zeros(eltype(a), k)
    i = 1
    while i <= k && length(h) > 0
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

"""
    CornerSpaceRenorm.corner_subspace(ρA, ρB, M)

Find the `SubspaceBasis` associated with the `M`-dimensional corner defined by the
pair of states `(ρA, ρB)`, the pairs of eigenkets of A and B defining the corner
basis and the eigenkets of A and B.
# Arguments
* `ρA`: density matrix of the subsystem A.
* `ρB`: density matrix of the subsystem B.
* `M`: corner dimension.
"""
function corner_subspace(ρA::DenseMatrix, bA::Basis, ρB::DenseMatrix, bB::Basis, M::Int)
    ps_A, αs_A = eigen(ρA); # αs_A[:,i] : i-th eigenvector
    ps_B, αs_B = eigen(ρB);
    mgs!(αs_A);
    mgs!(αs_B);
    @inbounds ϕs_A = [Ket(bA, αs_A[:,i]) for i in 1:length(ps_A)]
    @inbounds ϕs_B = [Ket(bB, αs_B[:,i]) for i in 1:length(ps_B)]
    handles, prod_pairs = max_prod_pairs(real.(ps_A), real.(ps_B), M)
    bC = typeof(bA) <: CompositeBasis ? CompositeBasis(bA.bases...,bB.bases...) : CompositeBasis(bA,bB)
    return CornerBasis(bA, bB, M), handles, ϕs_A, ϕs_B
end
function corner_subspace(ρA::AbstractOperator{B1,B1}, ρB::AbstractOperator{B2,B2}, M::Int) where {B1<:Basis,B2<:Basis}
    if isdense(ρA) && isdense(ρB)
        return corner_subspace(ρA.data, ρA.basis_l, ρB.data, ρB.basis_l, M)
    else
        return corner_subspace(Matrix(ρA.data), ρA.basis_l, Matrix(ρB.data), ρB.basis_l, M)
    end
end

"""
    hermitianize(x)

Extract the Hermitian part of `x`.
# Arguments
* `x`: some operator.
"""
function hermitianize(x::AbstractOperator{B,B}) where {B<:Basis}
    ishermitian(x) && return x
    n = LinearAlgebra.checksquare(x.data);
    y = copy(x);
    @inbounds @views for j = 1:n, i = 1:j
        y.data[i,j] = (y.data[i,j] + conj(y.data[j,i])) / 2.
        y.data[j,i] = conj(y.data[i,j]);
    end
    return y
end

"""
    hermitianize!(x)

Extract the Hermitian part of `x`.
# Arguments
* `x`: some operator.
"""
function hermitianize!(x::AbstractOperator{B,B}) where {B<:Basis}
    ishermitian(x) && return x
    n = LinearAlgebra.checksquare(x.data);
    @inbounds @views for j = 1:n, i = 1:j
        x.data[i,j] = (x.data[i,j] + conj(x.data[j,i])) / 2.
        x.data[j,i] = conj(x.data[i,j]);
    end
    return x
end

"""
    merge(s1, s2, d, ρ1, ρ2, M)

Merge two `NdSystem`s along the `d`-th direction with corner compression.
# Arguments
* `s1`: `NdSystem`.
* `s2`: `NdSystem`.
* `d`: merging direction.
* `ρ1`: state of the system `s1`.
* `ρ2`: state of the system `s2`.
* `M`: corner dimension.
"""
function Base.merge(s1::NdSystem{N},s2::NdSystem{N},d::Integer,
                    ρ1::Operator{B1,B1,T1},ρ2::Operator{B2,B2,T2},M::Int;
                    return_dm=false) where {N,B1<:Basis,B2<:Basis,T1<:DenseMatrix,T2<:DenseMatrix}
    # TO DO: tests
    @assert 0 < M <= length(s1.gbasis)*length(s2.gbasis) "The corner dimension must not exceed that of the total Hilbert space."
    length(s1.lattice.Vext[d]) == length(s2.lattice.Vint[d]) || throw(DimensionMismatch("Lattices cannot be merged in this direction."))
    @assert s1.lattice.pbc == s2.lattice.pbc "Cannot merge systems with periodic and open open boundary conditions."
    lattice = union(s1.lattice,s2.lattice,d)

    all(s1.trate .≈ s2.trate) || @info "The two input systems do not have equal tunnelling rates along all directions."

    bC, handles, ϕs_1, ϕs_2 = corner_subspace(ρ1,ρ2,M)
    vt_1 = map(x->transpose(x.data), ϕs_1)
    vc_1 = map(x->conj(x.data), ϕs_1)
    vt_2 = map(x->transpose(x.data), ϕs_2)
    vc_2 = map(x->conj(x.data), ϕs_2)
    cache1 = [similar(first(vc_1)) for i in 1:Threads.nthreads()]
    cache2 = [similar(first(vc_2)) for i in 1:Threads.nthreads()]

    function 𝒫1(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 2
        opC = DenseOperator(bC)
        Threads.@threads for j in 1:length(handles)
            hj = handles[j]
            mul!(cache1[Threads.threadid()], op.data, vc_1[hj[1]])
            for i in 1:length(handles)
                hi = handles[i]
                opC.data[i,j] = vt_1[hi[1]] * cache1[Threads.threadid()] * (vt_2[hi[2]]*vc_2[hj[2]])
            end
        end
        return opC
    end
    function 𝒫1_bis(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 2
        opC = DenseOperator(bC)
        for j in 1:length(handles)
            hj = handles[j]
            mul!(cache1[Threads.threadid()], op.data, vc_1[hj[1]])
            for i in 1:length(handles)
                hi = handles[i]
                opC.data[i,j] = vt_1[hi[1]] * cache1[Threads.threadid()] * (vt_2[hi[2]]*vc_2[hj[2]])
            end
        end
        return opC
    end
    function 𝒫2(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 1
        opC = DenseOperator(bC)
        Threads.@threads for j in 1:length(handles)
            hj = handles[j]
            mul!(cache2[Threads.threadid()], op.data, vc_2[hj[2]])
            for i in 1:length(handles)
                hi = handles[i]
                opC.data[i,j] = vt_2[hi[2]] * cache2[Threads.threadid()] * (vt_1[hi[1]]*vc_1[hj[1]])
            end
        end
        return opC
    end
    function 𝒫(op1,op2)
        opC = DenseOperator(bC)
        Threads.@threads for j in 1:length(handles)
            hj = handles[j]
            mul!(cache1[Threads.threadid()], op1.data, vc_1[hj[1]])
            mul!(cache2[Threads.threadid()], op2.data, vc_2[hj[2]])
            for i in 1:length(handles)
                hi = handles[i]
                opC.data[i,j] = (vt_1[hi[1]] * cache1[Threads.threadid()]) * (vt_2[hi[2]] * cache2[Threads.threadid()])
            end
        end
        return opC
    end

    # TO DO: exploit Hermicianity of H to compute half of the matrix elements in the corner
    H = DenseOperator(bC)
    if lattice.pbc
        Ht1 = zeros(ComplexF64,size(s1.H.data))
        if s1.lattice.shape[d] > 2
            @inbounds for i in 1:length(s1.Htext[d])
                Ht1 .+= s1.trate[d] * (s1.Htext[d][i] * dagger(s1.Htint[d][i])).data;
            end
        end
        Ht2 = zeros(ComplexF64,size(s2.H.data))
        if s2.lattice.shape[d] > 2
            @inbounds for i in 1:length(s2.Htext[d])
                Ht2 .+= s2.trate[d] * (s2.Htext[d][i] * dagger(s2.Htint[d][i])).data;
            end
        end
        H.data .= 𝒫1(DenseOperator(s1.gbasis, s1.H.data .- (Ht1 .+ Ht1'))).data .+ 𝒫2(DenseOperator(s2.gbasis, s2.H.data .- (Ht2 .+ Ht2'))).data;
    else
        H.data .= 𝒫1(s1.H).data .+ 𝒫2(s2.H).data;
    end
    gbasis = H.basis_l;
    Ht = zeros(ComplexF64,size(H.data));
    @inbounds for i in 1:length(s1.Htext[d])
        Ht .= s1.trate[d] * 𝒫(s1.Htext[d][i],dagger(s2.Htint[d][i])).data;
        if lattice.pbc
            Ht .+= s1.trate[d] * 𝒫(dagger(s1.Htint[d][i]), s2.Htext[d][i]).data;
        end
        H.data .+= Ht .+ Ht';
    end
    hermitianize!(H);

    J = [[𝒫1(s1.J[i]) for i in 1:length(s1.J)]; [𝒫2(s2.J[i]) for i in 1:length(s2.J)]]

    #return H
    d⊥ = [_d for _d in 1:N if _d!= d]
    T = typeof(H)
    Htint::Array{Array{T,1},1} = [Array{T,1}(undef,0) for i in 1:N]
    Htext::Array{Array{T,1},1} = [Array{T,1}(undef,0) for i in 1:N]

    Htint[d] = 𝒫1.(s1.Htint[d])
    Htext[d] = 𝒫2.(s2.Htext[d])
    for _d in d⊥
        Htint[_d] = [𝒫1.(s1.Htint[_d]); 𝒫2.(s2.Htint[_d])]
        Htext[_d] = [𝒫1.(s1.Htext[_d]); 𝒫2.(s2.Htext[_d])]
    end

    # TO DO: Test that observable types are compatible
    obs = [[Dict{String,typeof(H)}() for i in 1:nv(s1.lattice)]; [Dict{String,typeof(H)}() for i in 1:nv(s2.lattice)]]

    for i in 1:length(s1.observables)
        for k in keys(s1.observables[i])
            obs[i][k] = 𝒫1(s1.observables[i][k])
        end
    end
    for i in 1:length(s2.observables)
        for k in keys(s2.observables[i])
            obs[length(s1.observables)+i][k] = 𝒫2(s2.observables[i][k])
        end
    end

    if return_dm
        return NdSystem{N,typeof(lattice),typeof(gbasis),typeof(H),eltype(first(Htint)),eltype(J),typeof(H)}(lattice,gbasis,H,s1.trate,Tuple(Htint),Tuple(Htext),J,obs), normalize(𝒫(ρ1, ρ2))
    else
        return NdSystem{N,typeof(lattice),typeof(gbasis),typeof(H),eltype(first(Htint)),eltype(J),typeof(H)}(lattice,gbasis,H,s1.trate,Tuple(Htint),Tuple(Htext),J,obs)
    end
end
