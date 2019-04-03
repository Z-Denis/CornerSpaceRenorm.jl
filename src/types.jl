using LightGraphs

const AbstractLattice = AbstractGraph{Int64}
abstract type Lattice <: AbstractLattice end
abstract type AbstractCorner end
abstract type AbstractSystem end

struct SquareLattice <: Lattice
    L::SimpleGraph{Int64}
    nx::Int64
    ny::Int64

    extsites::Vector{Int64}
end

function SquareLattice(nx::Integer,ny::Integer,lbasis::B;periodic::Bool = false) where B<:Basis
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

struct System{B<:Basis,L<:Lattice,N,
              O1<:AbstractOperator{B,B},
              O2<:AbstractOperator{B,B},
              O3<:AbstractOperator{B,B}} <: AbstractSystem
    # Lattice
    lattice::L
    # Bases
    lbasis::B
    gbasis::CompositeBasis{Tuple{Vararg{FockBasis,N}}} # N == nx * ny
    # Operators
    lH::O1
    tH::Tuple{O2,O2}
    J::Vector{O3}
end

function System(lat::L,lH::O1,tH::Tuple{O2,O2},J::Vector{O3}) where {B<:Basis,L<:Lattice,
                                                                     O1<:AbstractOperator{B,B},
                                                                     O2<:AbstractOperator{B,B},
                                                                     O3<:AbstractOperator{B,B}}
    @assert lH.basis_r == lH.basis_l
    @assert all([tH[i].basis_l == lH.basis_l for i in 1:2])
    @assert all([tH[i].basis_r == lH.basis_l for i in 1:2])
    @assert all([J[i].basis_l == lH.basis_l for i in 1:length(J)])
    @assert all([J[i].basis_r == lH.basis_l for i in 1:length(J)])
    return System{B,L,nv(lat),O1,O2,O3}(lat,lH.basis_l,CompositeBasis([lH.basis_l for i in 1:nv(lat)]...),lH,tH,J)
end
