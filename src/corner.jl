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
    return System{typeof(s.lattice),typeof(cspace),typeof(H),eltype(Httop),eltype(J)}(s.lattice,cspace,H,Httop,Htbottom,Htleft,Htright,J)
end

function merge(s1::AbstractSystem,s2::AbstractSystem)
    if s1.lattice.ny == s2.lattice.ny
        return vmerge(s1,s2)
    elseif s1.lattice.ny == s2.lattice.ny
        return hmerge(s1,s2)
    else
        error("Systems have incompatible lattice sizes: ($(s1.lattice.nx),$(s1.lattice.ny)) and ($(s2.lattice.nx),$(s2.lattice.ny))")
    end
end

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

    return System{typeof(lattice),typeof(gbasis),typeof(H),eltype(Httop),eltype(J)}(lattice,gbasis,H,Httop,Htbottom,Htleft,Htright,J)
end

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

    return System{typeof(lattice),typeof(gbasis),typeof(H),eltype(Httop),eltype(J)}(lattice,gbasis,H,Httop,Htbottom,Htleft,Htright,J)
end
