using DataStructures, LinearAlgebra

"""
    maxk(a, k)

Find `k` largest elements of an array and return their indices.
# Arguments
* `a`: some array.
* `k`: number of largest elements to be found.
"""
function maxk(a::AbstractVector, k::Int)
    ix = partialsortperm(a, 1:k, rev=true)
    @inbounds @views collect(ix), collect(a[ix[1:k]])
end

"""
    max_prod_pairs(a, b, k)

Find `k` pairs `(a_i,b_i)` with largest products from two arrays `a` and `b`.
# Arguments
* `a`: some array.
* `b`: some array.
* `k`: number of largest products to be found.
"""
function max_prod_pairs(a::AbstractVector, b::AbstractVector, k::Int)
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

"""
    corner_subspace(ρA, ρB, M)

Find the `SubspaceBasis` associated with the `M`-dimensional corner defined by the
pair of states `(ρA, ρB)`, the pairs of eigenkets of A and B defining the corner
basis and the eigenkets of A and B.
# Arguments
* `ρA`: density matrix of the subsystem A.
* `ρB`: density matrix of the subsystem B.
* `M`: corner dimension.
"""
function corner_subspace(ρA::DenseOperator{B1,B1}, ρB::DenseOperator{B2,B2}, M::Int) where {B1<:Basis,B2<:Basis}
    bA = ρA.basis_l
    bB = ρB.basis_l
    ps_A, αs_A = eigen(ρA.data); # αs_A[:,i] : i-th eigenvector
    ps_B, αs_B = eigen(ρB.data);
    @inbounds ϕs_A = [Ket(bA, αs_A[:,i]) for i in 1:length(ps_A)]
    @inbounds ϕs_B = [Ket(bB, αs_B[:,i]) for i in 1:length(ps_B)]
    handles, prod_pairs = max_prod_pairs(real.(ps_A), real.(ps_B), M)
    bC = typeof(bA) <: CompositeBasis ? CompositeBasis(bA.bases...,bB.bases...) : CompositeBasis(bA,bB)
    return SubspaceBasis(bC, [ϕs_A[idcs[1]] ⊗ ϕs_B[idcs[2]] for idcs in handles]), handles, ϕs_A, ϕs_B
end

"""
    cornerize(s, cspace)

Project a system into a given corner space.
# Arguments
* `s`: `System`.
* `cspace`: `SubspaceBasis` of the corner.
"""
function cornerize(s::AbstractSystem,cspace::SubspaceBasis)
    P = projector(cspace,s.gbasis);
    Pd = dagger(P);
    proj(op) = P*op*Pd;
    H = proj(s.H);
    Httop = proj.(s.Httop);
    Htbottom = proj.(s.Htbottom);
    Htleft = proj.(s.Htleft);
    Htright = proj.(s.Htright);
    J = proj.(s.J);
    obs = [Dict([name => proj(lop) for (name,lop) in s.observables[i]]) for i in 1:nv(s.lattice)]
    return System{typeof(s.lattice),typeof(cspace),typeof(H),eltype(Httop),eltype(J),typeof(H)}(s.lattice,cspace,H,Httop,Htbottom,Htleft,Htright,J,obs)
end

"""
    merge(s1, s2)

Merge two `System`s along some compatible dimension with no corner compression.
# Arguments
* `s1`: `System`.
* `s2`: `System`.
"""
function merge(s1::AbstractSystem,s2::AbstractSystem)
    if s1.lattice.ny == s2.lattice.ny
        return vmerge(s1,s2)
    elseif s1.lattice.ny == s2.lattice.ny
        return hmerge(s1,s2)
    else
        error("Systems have incompatible lattice sizes: ($(s1.lattice.nx),$(s1.lattice.ny)) and ($(s2.lattice.nx),$(s2.lattice.ny))")
    end
end

"""
    vmerge(s1, s2)

Merge two `System`s vertically with no corner compression.
# Arguments
* `s1`: `System`.
* `s2`: `System`.
"""
function vmerge(s1::AbstractSystem,s2::AbstractSystem)
    # TO DO: tests
    lattice = vunion(s1.lattice,s2.lattice)
    Id1 = one(s1.gbasis)
    Id2 = one(s2.gbasis)
    H = s1.H ⊗ Id2 + Id1 ⊗ s2.H
    gbasis = H.basis_l;

    @inbounds for i in 1:length(s1.Htbottom)
        H.data .+= (s1.Htbottom[i] ⊗ dagger(s2.Httop[i])).data;
        H.data .+= (dagger(s1.Htbottom[i]) ⊗ s2.Httop[i]).data;
    end
    J = ([s1.J[i] ⊗ Id2 for i in 1:length(s1.J)]) ∪ ([Id1 ⊗ s2.J[i] for i in 1:length(s2.J)])
    Httop = [s1.Httop[i] ⊗ Id2 for i in 1:length(s1.Httop)]
    Htbottom = [Id1 ⊗ s2.Htbottom[i] for i in 1:length(s2.Htbottom)]
    Htleft = [s1.Htleft[i] ⊗ Id2 for i in 1:length(s1.Htleft)] ∪ [Id1 ⊗ s2.Htleft[i] for i in 1:length(s2.Htleft)]
    Htright = [s1.Htright[i] ⊗ Id2 for i in 1:length(s1.Htright)] ∪ [Id1 ⊗ s2.Htright[i] for i in 1:length(s2.Htright)]

    # TO DO: Test that observable types are compatible
    obs = [[Dict{String,typeof(H)}() for i in 1:nv(s1.lattice)]; [Dict{String,typeof(H)}() for i in 1:nv(s2.lattice)]]

    for i in 1:length(s1.observables)
        for k in keys(s1.observables[i])
            obs[i][k] = s1.observables[i][k] ⊗ Id2
        end
    end
    for i in 1:length(s2.observables)
        for k in keys(s2.observables[i])
            obs[length(s1.observables)+i][k] = Id1 ⊗ s2.observables[i][k]
        end
    end
    return System{typeof(lattice),typeof(gbasis),typeof(H),eltype(Httop),eltype(J),typeof(H)}(lattice,gbasis,H,Httop,Htbottom,Htleft,Htright,J,obs)
end

"""
    vmerge(s1, s2)

Merge two `System`s horizontally with no corner compression.
# Arguments
* `s1`: `System`.
* `s2`: `System`.
"""
function hmerge(s1::AbstractSystem,s2::AbstractSystem)
    # TO DO: tests
    lattice = hunion(s1.lattice,s2.lattice)
    Id1 = one(s1.gbasis)
    Id2 = one(s2.gbasis)
    H = s1.H ⊗ Id2 + Id1 ⊗ s2.H
    gbasis = H.basis_l;

    @inbounds for i in 1:length(s1.Htright)
        H.data .+= (s1.Htright[i] ⊗ dagger(s2.Htleft[i])).data;
        H.data .+= (dagger(s1.Htright[i]) ⊗ s2.Htleft[i]).data;
    end
    J = ([s1.J[i] ⊗ Id2 for i in 1:length(s1.J)]) ∪ ([Id1 ⊗ s2.J[i] for i in 1:length(s2.J)])
    Httop = [s1.Httop[i] ⊗ Id2 for i in 1:length(s1.Httop)] ∪ [Id1 ⊗ s2.Httop[i] for i in 1:length(s2.Httop)]
    Htbottom = [s1.Htbottom[i] ⊗ Id2 for i in 1:length(s1.Htbottom)] ∪ [Id1 ⊗ s2.Htbottom[i] for i in 1:length(s2.Htbottom)]
    Htleft = [s1.Htleft[i] ⊗ Id2 for i in 1:length(s1.Htleft)]
    Htright = [Id1 ⊗ s2.Htright[i] for i in 1:length(s2.Htright)]

    # TO DO: Test that observable types are compatible
    obs = [[Dict{String,typeof(H)}() for i in 1:nv(s1.lattice)]; [Dict{String,typeof(H)}() for i in 1:nv(s2.lattice)]]

    for i in 1:length(s1.observables)
        for k in keys(s1.observables[i])
            obs[i][k] = s1.observables[i][k] ⊗ Id2
        end
    end
    for i in 1:length(s2.observables)
        for k in keys(s2.observables[i])
            obs[length(s1.observables)+i][k] = Id1 ⊗ s2.observables[i][k]
        end
    end
    return System{typeof(lattice),typeof(gbasis),typeof(H),eltype(Httop),eltype(J),typeof(H)}(lattice,gbasis,H,Httop,Htbottom,Htleft,Htright,J,obs)
end

"""
    hermitianize!(x)

Extract the Hermitian part of `x`.
# Arguments
* `x`: some operator.
"""
function hermitianize!(x::AbstractOperator{B,B}) where {B<:Basis}
    n = LinearAlgebra.checksquare(x.data);
    @inbounds for j = 2:n, i = 1:j-1
        x.data[i,j] = (x.data[i,j] + conj(x.data[j,i])) / 2
        x.data[j,i] = conj(x.data[i,j]);
    end
    nothing
end

"""
    merge(s1, s2, ρ1, ρ2, M)

Merge two `System`s along some compatible dimension with corner compression.
# Arguments
* `s1`: `System`.
* `s2`: `System`.
* `ρ1`: state of the system `s1`.
* `ρ2`: state of the system `s2`.
* `M`: corner dimension.
"""
function Base.merge(s1::AbstractSystem,s2::AbstractSystem,ρ1::DenseOperator{B1,B1},ρ2::DenseOperator{B2,B2},M::Int) where {B1<:Basis,B2<:Basis}
    # TO DO: add tests on M
    if s1.lattice.ny == s2.lattice.ny
        return vmerge(s1,s2,ρ1,ρ2,M)
    elseif s1.lattice.ny == s2.lattice.ny
        return hmerge(s1,s2,ρ1,ρ2,M)
    else
        error("Systems have incompatible lattice sizes: ($(s1.lattice.nx),$(s1.lattice.ny)) and ($(s2.lattice.nx),$(s2.lattice.ny))")
    end
end

"""
    vmerge(s1, s2, ρ1, ρ2, M)

Merge two `System`s vertically along some compatible dimension with corner compression.
# Arguments
* `s1`: `System`.
* `s2`: `System`.
* `ρ1`: state of the system `s1`.
* `ρ2`: state of the system `s2`.
* `M`: corner dimension.
"""
function vmerge(s1::AbstractSystem,s2::AbstractSystem,ρ1::DenseOperator{B1,B1},ρ2::DenseOperator{B2,B2},M::Int) where {B1<:Basis,B2<:Basis}
    # TO DO: tests
    lattice = vunion(s1.lattice,s2.lattice)

    bC, handles, ϕs_1, ϕs_2 = corner_subspace(ρ1,ρ2,M)
    function 𝒫1(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 2
        return DenseOperator(bC,[transpose(ϕs_1[hi[1]].data) * (op * conj.(ϕs_1[hj[1]])).data * (transpose(ϕs_2[hi[2]].data) * conj.(ϕs_2[hj[2]]).data) for hi in handles, hj in handles])
    end
    function 𝒫2(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 1
        return DenseOperator(bC,[transpose(ϕs_2[hi[2]].data) * (op * conj.(ϕs_2[hj[2]])).data * (transpose(ϕs_1[hi[1]].data) * conj.(ϕs_1[hj[1]]).data) for hi in handles, hj in handles])
    end
    function 𝒫(op1,op2)
        return DenseOperator(bC,[(transpose(ϕs_1[hi[1]].data) * (op1 * conj.(ϕs_1[hj[1]])).data) * (transpose(ϕs_2[hi[2]].data) * (op2 * conj.(ϕs_2[hj[2]])).data) for hi in handles, hj in handles])
    end

    # TO DO: exploit Hermicianity of H to compute half of the matrix elements in the corner
    H = 𝒫1(s1.H) + 𝒫2(s2.H);
    gbasis = H.basis_l;
    @inbounds for i in 1:length(s1.Htbottom)
        Ht = 𝒫(s1.Htbottom[i],dagger(s2.Httop[i])).data;
        H.data .+= Ht .+ Ht';
    end
    hermitianize!(H);

    J = [𝒫1(s1.J[i]) for i in 1:length(s1.J)] ∪ [𝒫2(s2.J[i]) for i in 1:length(s2.J)]

    Httop = [𝒫1(s1.Httop[i]) for i in 1:length(s1.Httop)];
    Htbottom = [𝒫2(s2.Htbottom[i]) for i in 1:length(s2.Htbottom)];
    Htleft = [𝒫1(s1.Htleft[i]) for i in 1:length(s1.Htleft)] ∪ [𝒫2(s2.Htleft[i]) for i in 1:length(s2.Htleft)];
    Htright = [𝒫1(s1.Htright[i]) for i in 1:length(s1.Htright)] ∪ [𝒫2(s2.Htright[i]) for i in 1:length(s2.Htright)];

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

    return System{typeof(lattice),typeof(gbasis),typeof(H),eltype(Httop),eltype(J),typeof(H)}(lattice,gbasis,H,Httop,Htbottom,Htleft,Htright,J,obs)
end

"""
    hmerge(s1, s2, ρ1, ρ2, M)

Merge two `System`s horizontally along some compatible dimension with corner compression.
# Arguments
* `s1`: `System`.
* `s2`: `System`.
* `ρ1`: state of the system `s1`.
* `ρ2`: state of the system `s2`.
* `M`: corner dimension.
"""
function hmerge(s1::AbstractSystem,s2::AbstractSystem,ρ1::DenseOperator{B1,B1},ρ2::DenseOperator{B2,B2},M::Int) where {B1<:Basis,B2<:Basis}
    # TO DO: tests
    lattice = hunion(s1.lattice,s2.lattice)

    bC, handles, ϕs_1, ϕs_2 = corner_subspace(ρ1,ρ2,M)
    function 𝒫1(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 2
        return DenseOperator(bC,[transpose(ϕs_1[hi[1]].data) * (op * conj.(ϕs_1[hj[1]])).data * (transpose(ϕs_2[hi[2]].data) * conj.(ϕs_2[hj[2]]).data) for hi in handles, hj in handles])
    end
    function 𝒫2(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 1
        return DenseOperator(bC,[transpose(ϕs_2[hi[2]].data) * (op * conj.(ϕs_2[hj[2]])).data * (transpose(ϕs_1[hi[1]].data) * conj.(ϕs_1[hj[1]]).data) for hi in handles, hj in handles])
    end
    function 𝒫(op1,op2)
        return DenseOperator(bC,[(transpose(ϕs_1[hi[1]].data) * (op1 * conj.(ϕs_1[hj[1]])).data) * (transpose(ϕs_2[hi[2]].data) * (op2 * conj.(ϕs_2[hj[2]])).data) for hi in handles, hj in handles])
    end

    # TO DO: exploit Hermicianity of H to compute half of the matrix elements in the corner
    H = 𝒫1(s1.H) + 𝒫2(s2.H);
    gbasis = H.basis_l;
    @inbounds for i in 1:length(s1.Htright)
        Ht = 𝒫(s1.Htright[i],dagger(s2.Htleft[i])).data;
        H.data .+= Ht .+ Ht';
    end
    hermitianize!(H);
    J = [𝒫1(s1.J[i]) for i in 1:length(s1.J)] ∪ [𝒫2(s2.J[i]) for i in 1:length(s2.J)]

    Httop = [𝒫1(s1.Httop[i]) for i in 1:length(s1.Httop)] ∪ [𝒫2(s2.Httop[i]) for i in 1:length(s2.Httop)];
    Htbottom = [𝒫1(s1.Htbottom[i]) for i in 1:length(s1.Htbottom)] ∪ [𝒫2(s2.Htbottom[i]) for i in 1:length(s2.Htbottom)];
    Htleft = [𝒫1(s1.Htleft[i]) for i in 1:length(s1.Htleft)];
    Htright = [𝒫2(s2.Htright[i]) for i in 1:length(s2.Htright)];

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

    return System{typeof(lattice),typeof(gbasis),typeof(H),eltype(Httop),eltype(J),typeof(H)}(lattice,gbasis,H,Httop,Htbottom,Htleft,Htright,J,obs)
end

"""
    merge(s1, s2, d, ρ1, ρ2, M)

Merge two `ZnSystem`s along the `d`-th direction with corner compression.
# Arguments
* `s1`: `ZnSystem`.
* `s2`: `ZnSystem`.
* `d`: merging direction.
* `ρ1`: state of the system `s1`.
* `ρ2`: state of the system `s2`.
* `M`: corner dimension.
"""
function Base.merge(s1::ZnSystem{N},s2::ZnSystem{N},d::Integer,ρ1::DenseOperator{B1,B1},ρ2::DenseOperator{B2,B2},M::Int) where {N,B1<:Basis,B2<:Basis}
    # TO DO: tests
    lattice = union(s1.lattice,s2.lattice,d)

    bC, handles, ϕs_1, ϕs_2 = corner_subspace(ρ1,ρ2,M)
    function 𝒫1(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 2
        return DenseOperator(bC,[transpose(ϕs_1[hi[1]].data) * (op * conj.(ϕs_1[hj[1]])).data * (transpose(ϕs_2[hi[2]].data) * conj.(ϕs_2[hj[2]]).data) for hi in handles, hj in handles])
    end
    function 𝒫2(op)
        # TO DO: take advantage of orthogonormality to get rid of the scalar product on subspace 1
        return DenseOperator(bC,[transpose(ϕs_2[hi[2]].data) * (op * conj.(ϕs_2[hj[2]])).data * (transpose(ϕs_1[hi[1]].data) * conj.(ϕs_1[hj[1]]).data) for hi in handles, hj in handles])
    end
    function 𝒫(op1,op2)
        return DenseOperator(bC,[(transpose(ϕs_1[hi[1]].data) * (op1 * conj.(ϕs_1[hj[1]])).data) * (transpose(ϕs_2[hi[2]].data) * (op2 * conj.(ϕs_2[hj[2]])).data) for hi in handles, hj in handles])
    end

    # TO DO: exploit Hermicianity of H to compute half of the matrix elements in the corner
    H = 𝒫1(s1.H) + 𝒫2(s2.H);
    gbasis = H.basis_l;
    @inbounds for i in 1:length(s1.Htext[d])
        Ht = 𝒫(s1.Htext[d][i],dagger(s2.Htint[d][i])).data;
        H.data .+= Ht .+ Ht';
    end
    hermitianize!(H);

    J = [𝒫1(s1.J[i]) for i in 1:length(s1.J)] ∪ [𝒫2(s2.J[i]) for i in 1:length(s2.J)]

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

    return ZnSystem{N,typeof(lattice),typeof(gbasis),typeof(H),eltype(first(Htint)),eltype(J),typeof(H)}(lattice,gbasis,H,Tuple(Htint),Tuple(Htext),J,obs)
end
