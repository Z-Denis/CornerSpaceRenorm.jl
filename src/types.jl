using LightGraphs

const AbstractLattice = AbstractGraph{Int64}
abstract type Lattice <: AbstractLattice end

struct SquareLattice <: Lattice
    L::SimpleGraph{Int64}
    nx::Int64
    ny::Int64
    
    extsites::Vector{Int64}
end

function SquareLattice(nx::Integer,ny::Integer;periodic::Bool = false)
    L = Grid([nx, ny];periodic=periodic);
    V = vertices(L)
    extsites = [V[1:nx]; V[end-nx+1:end];[V[nx*i] for i in 2:ny-1];[V[nx*i+1] for i in 1:ny-2]];
    return SquareLattice(L,nx,ny,extsites)
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
@inline extsites(L::Lattice) = L.extsites;

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
    tH::Tuple{O2,O2}
    J::Vector{O3}
end

function System(lat::L,H::O1,tH::Tuple{O2,O2},J::Vector{O3}) where {L<:Lattice,B<:Basis,
                                                                    O1<:AbstractOperator{B,B},
                                                                    O2<:AbstractOperator{B,B},
                                                                    O3<:AbstractOperator{B,B}}
    gbasis = H.basis_l;
    @assert H.basis_r == H.basis_l == gbasis
    @assert first(tH).basis_l == gbasis
    @assert first(tH).basis_r == gbasis
    @assert first(J).basis_l == gbasis
    @assert first(J).basis_r == gbasis
    return System{L,B,O1,O2,O3}(lat,gbasis,H,tH,J)
end

struct CornerBasis{T<:Ket} <: Basis
    shape::Vector{Int}
    M::Int
    basisstates::Vector{T}
    function CornerBasis(kets::Vector{K}) where K<:Ket
        M = length(kets);
        new{eltype(kets)}([M], M, kets);
    end
end
