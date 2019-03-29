using LightGraphs

const AbstractLattice = AbstractGraph{Int64}
abstract type Lattice <: AbstractLattice end
abstract type AbstractCorner end

struct SquareLattice{B<:Basis,N} <: Lattice
    L::SimpleGraph{Int64}
    nx::Int64
    ny::Int64
    lbasis::B
    gbasis::CompositeBasis{Tuple{Vararg{FockBasis,N}}} # N == nx * ny
end

function SquareLattice(nx::Integer,ny::Integer,lbasis::B;periodic::Bool = false) where B<:Basis
    L = Grid([nx, ny];periodic=periodic);
    gbasis = CompositeBasis([lbasis for i in 1:LightGraphs.nv(L)]...)
    return SquareLattice{B,LightGraphs.nv(L)}(L,nx,ny,lbasis,gbasis)
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

#=
function vunion(L1::SquareLattice,L2::SquareLattice)
    L1_bottom = [L1.nx*i for i in 1:L1.ny]
    L2_top = [L2.nx*(i-1)+1 for i in 1:L2.ny]
    L = blockdiag(L1.L,L2.L);
    for i in L1.ny
        add_edge!(L,Edge(L1_bottom[i],nv(L1)+L2_top[i])
    end
    return SquareLattice(L1.nx+L2.nx,L1.ny,L1.lbasis)
end
=#
function Base.:union(L1::SquareLattice,L2::SquareLattice)
    @assert L1.lbasis == L2.lbasis
    if L1.ny == L2.ny
        return vunion(L1,L2)
    elseif L1.nx == L2.nx
        return hunion(L1,L2)
    else
        error("Incompatible lattice sizes: ($(L1.nx),$(L1.ny)) and ($(L2.nx),$(L2.ny))")
    end
end

function vunion(L1::SquareLattice,L2::SquareLattice)
    return SquareLattice(L1.nx+L2.nx,L1.ny,L1.lbasis)
end

function hunion(L1::SquareLattice,L2::SquareLattice)
    return SquareLattice(L1.nx,L1.ny+L2.ny,L1.lbasis)
end
