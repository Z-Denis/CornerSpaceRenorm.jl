using LightGraphs

const AbstractLattice = AbstractGraph{Int64}
abstract type Lattice <: AbstractLattice end

struct SquareLattice <: Lattice
    L::SimpleGraph{Int64}
    nx::Int64
    ny::Int64

    Vtop::Vector{Int64}
    Vbottom::Vector{Int64}
    Vleft::Vector{Int64}
    Vright::Vector{Int64}
end

function SquareLattice(nx::Integer,ny::Integer;periodic::Bool = false)
    L = Grid([nx, ny];periodic=periodic);
    V = vertices(L)
    Vtop    = [V[nx*i+1] for i in 0:ny-1];
    Vbottom = [V[nx*i] for i in 1:ny];
    Vleft   = V[1:nx];
    Vright  = V[end-nx+1:end];
    return SquareLattice(L,nx,ny,Vtop,Vbottom,Vleft,Vright)
end

@inline Base.eltype(L::Lattice) = Base.eltype(L.L);
@inline LightGraphs.vertices(L::Lattice) = LightGraphs.vertices(L.L);
@inline LightGraphs.edges(L::Lattice) = LightGraphs.edges(L.L);
@inline LightGraphs.nv(L::Lattice) = LightGraphs.nv(L.L);
@inline LightGraphs.ne(L::Lattice) = LightGraphs.ne(L.L);
@inline LightGraphs.has_edge(L::Lattice) = LightGraphs.has_edge(L.L);
@inline LightGraphs.has_vertex(L::Lattice) = LightGraphs.has_vertex(L.L);
@inline LightGraphs.inneighbors(L::Lattice) = LightGraphs.inneighbors(L.L);
@inline LightGraphs.outneighbors(L::Lattice) = LightGraphs.outneighbors(L.L);
@inline LightGraphs.is_directed(L::Lattice) = false;
@inline LightGraphs.is_directed(T::Type{L}) where {L<:Lattice} = false;
@inline extsites(L::Lattice) = [L.Vu; L.Vb; L.Vl[2:end-1]; L.Vr[2:end-1]];

function Base.:union(L1::SquareLattice,L2::SquareLattice)
    if L1.ny == L2.ny
        return vunion(L1,L2)
    elseif L1.nx == L2.nx
        return hunion(L1,L2)
    else
        error("Incompatible lattice sizes: ($(L1.nx),$(L1.ny)) and ($(L2.nx),$(L2.ny))")
    end
end

function vunion(L1::Lattice,L2::Lattice)
    L = blockdiag(L1.L,L2.L);
    for i in 1:L1.ny
        add_edge!(L,Edge(L1.Vbottom[i],nv(L1)+L2.Vtop[i]))
    end

    return SquareLattice(L,L1.nx+L2.nx,L1.ny,L1.Vtop, nv(L1).+L2.Vbottom,
                                             L1.Vleft ∪ (nv(L1).+L2.Vleft),
                                             L1.Vright ∪ (nv(L1).+L2.Vright))
end

function hunion(L1::SquareLattice,L2::SquareLattice)
    L = blockdiag(L1.L,L2.L);
    for i in 1:L1.nx
        add_edge!(L,Edge(L1.Vright[i],nv(L1)+L2.Vleft[i]))
    end

    return SquareLattice(L,L1.nx,L1.ny+L2.ny,L1.Vtop ∪ (nv(L1).+L2.Vtop),
                                             L1.Vbottom ∪ (nv(L1).+L2.Vbottom),
                                             L1.Vleft, nv(L1).+L2.Vright)
end

abstract type AbstractSystem end

struct System{L<:Lattice,B<:Basis,
              O1<:AbstractOperator{B,B},O2<:AbstractOperator{B,B},
              O3<:AbstractOperator{B,B}} <: AbstractSystem
    # Lattice
    lattice::L
    # Bases
    gbasis::B
    # Operators
    H::O1
    Httop::Vector{O2}
    Htbottom::Vector{O2}
    Htleft::Vector{O2}
    Htright::Vector{O2}
    J::Vector{O3}
end

function System(lat::L,H::O1,lHt::O2,J::Vector{O3}) where {N,L<:Lattice,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B}}
    gbasis = H.basis_l;
    @assert nv(lat) == length(gbasis.bases)
    @assert lHt.basis_l == first(gbasis.bases)
    #@assert length(Htint) == 2lat.nx+2lat.ny-4
    #@assert length(Htext) == 2lat.nx+2lat.ny-4
    @assert H.basis_r == H.basis_l == gbasis
    #@assert first(Htint).basis_l == gbasis
    #@assert first(Htext).basis_r == gbasis
    @assert first(J).basis_l == gbasis
    @assert first(J).basis_r == gbasis
    Httop    = [embed(gbasis,[i],[lHt]) for i in lat.Vtop]
    Htbottom = [embed(gbasis,[i],[lHt]) for i in lat.Vbottom]
    Htleft   = [embed(gbasis,[i],[lHt]) for i in lat.Vleft]
    Htright  = [embed(gbasis,[i],[lHt]) for i in lat.Vright]

    return System{L,B,O1,eltype(Httop),O3}(lat,gbasis,H,Httop,Htbottom,Htleft,Htright,J)
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
    for i in 1:length(s1.Htbottom)
        Ht = s1.Htbottom[i] ⊗ dagger(s2.Httop[i]);
        H += Ht + dagger(Ht);
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
    for i in 1:length(s1.Htright)
        Ht = s1.Htright[i] ⊗ dagger(s2.Htleft[i]);
        H += Ht + dagger(Ht);
    end
    J = ([s1.J[i] ⊗ Id2 for i in 1:length(s1.J)]) ∪ ([Id1 ⊗ s2.J[i] for i in 1:length(s2.J)])
    Httop = [s1.Httop[i] ⊗ Id2 for i in 1:length(s1.Httop)] ∪ [Id1 ⊗ s2.Httop[i] for i in 1:length(s2.Httop)]
    Htbottom = [s1.Htbottom[i] ⊗ Id2 for i in 1:length(s1.Htbottom)] ∪ [Id1 ⊗ s2.Htbottom[i] for i in 1:length(s2.Htbottom)]
    Htleft = [s1.Htleft[i] ⊗ Id2 for i in 1:length(s1.Htleft)]
    Htright = [Id1 ⊗ s2.Htright[i] for i in 1:length(s2.Htright)]

    return System{typeof(lattice),typeof(gbasis),typeof(H),eltype(Httop),eltype(J)}(lattice,gbasis,H,Httop,Htbottom,Htleft,Htright,J)
end
