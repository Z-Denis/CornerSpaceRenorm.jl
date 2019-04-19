using LightGraphs

const AbstractLattice = AbstractGraph{Int64}
"""
    Lattice <: AbstractLattice

Abstract supertype for all lattices.
"""
abstract type Lattice <: AbstractLattice end

"""
    SquareLattice <: Lattice

Square lattice type.
"""
struct SquareLattice <: Lattice
    L::SimpleGraph{Int64}
    nx::Int64
    ny::Int64

    Vtop::Vector{Int64}
    Vbottom::Vector{Int64}
    Vleft::Vector{Int64}
    Vright::Vector{Int64}
end

"""
    SquareLattice(nx, ny; periodic = false)

Construct a `SquareLattice` of size `nx`*`ny`.

# Arguments
* `nx`: size of the first dimension.
* `ny`: size of the second dimension.
* `periodic=false`: boundary conditions.
"""
function SquareLattice(nx::Integer,ny::Integer;periodic::Bool = false)
    L = Grid([nx, ny];periodic=periodic);
    V = vertices(L)
    Vtop    = [V[nx*i+1] for i in 0:ny-1];
    Vbottom = [V[nx*i] for i in 1:ny];
    Vleft   = V[1:nx];
    Vright  = V[end-nx+1:end];
    return SquareLattice(L,nx,ny,Vtop,Vbottom,Vleft,Vright)
end

# Extend AbstractGraph{Int64} methods to all Lattice subtypes.
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

"""
    ∪(L1, L2)

Merge two `SquareLattice`s along some compatible dimension.
"""
function Base.:union(L1::SquareLattice,L2::SquareLattice)
    if L1.ny == L2.ny
        return vunion(L1,L2)
    elseif L1.nx == L2.nx
        return hunion(L1,L2)
    else
        error("Incompatible lattice sizes: ($(L1.nx),$(L1.ny)) and ($(L2.nx),$(L2.ny))")
    end
end

"""
    vunion(L1, L2)

Merge two `SquareLattice`s `L1` and `L2` vertically.
"""
function vunion(L1::Lattice,L2::Lattice)
    L = blockdiag(L1.L,L2.L);
    for i in 1:L1.ny
        add_edge!(L,Edge(L1.Vbottom[i],nv(L1)+L2.Vtop[i]))
    end

    return SquareLattice(L,L1.nx+L2.nx,L1.ny,L1.Vtop, nv(L1).+L2.Vbottom,
                                             L1.Vleft ∪ (nv(L1).+L2.Vleft),
                                             L1.Vright ∪ (nv(L1).+L2.Vright))
end

"""
    hunion(L1, L2)

Merge two `SquareLattice`s `L1` and `L2` horizontally.
"""
function hunion(L1::SquareLattice,L2::SquareLattice)
    L = blockdiag(L1.L,L2.L);
    for i in 1:L1.nx
        add_edge!(L,Edge(L1.Vright[i],nv(L1)+L2.Vleft[i]))
    end

    return SquareLattice(L,L1.nx,L1.ny+L2.ny,L1.Vtop ∪ (nv(L1).+L2.Vtop),
                                             L1.Vbottom ∪ (nv(L1).+L2.Vbottom),
                                             L1.Vleft, nv(L1).+L2.Vright)
end

"""
    AbstractSystem

Abstract supertype for all systems.
"""
abstract type AbstractSystem end

"""
    System <: AbstractSystem

Quantum dissipative system defined on a square lattice.
"""
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

"""
    System(lat, H, lHt, J)

Construct a `System` from a `Lattice` and operators defining the model.
# Arguments
* `lat`: Lattice.
* `H`: Hamiltonian of the full system.
* `lHt`: local tunneling operator.
* `J`: jump operators of the full system.
"""
function System(lat::L,H::O1,lHt::O2,J::Vector{O3}) where {N,L<:Lattice,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B}}
    gbasis = H.basis_l;
    @assert nv(lat) == length(gbasis.bases)
    @assert lHt.basis_l == first(gbasis.bases)
    @assert H.basis_r == H.basis_l == gbasis
    @assert first(J).basis_l == gbasis
    @assert first(J).basis_r == gbasis
    Httop    = [embed(gbasis,[i],[lHt]) for i in lat.Vtop]
    Htbottom = [embed(gbasis,[i],[lHt]) for i in lat.Vbottom]
    Htleft   = [embed(gbasis,[i],[lHt]) for i in lat.Vleft]
    Htright  = [embed(gbasis,[i],[lHt]) for i in lat.Vright]

    return System{L,B,O1,eltype(Httop),O3}(lat,gbasis,H,Httop,Htbottom,Htleft,Htright,J)
end
