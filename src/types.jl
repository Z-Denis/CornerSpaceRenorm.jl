using LightGraphs, GraphPlot

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

"""
    ZnLattice{N} <: Lattice

Z^`N` lattice type. Generic lattice with `N` directions.
"""
struct ZnLattice{N} <: Lattice
    L::SimpleGraph{Int64}
    shape::Tuple{Vararg{Int,N}}

    Vint::Tuple{Vararg{Vector{Int64},N}}
    Vext::Tuple{Vararg{Vector{Int64},N}}

    # Boundary conditions
    pbc::Bool
end

"""
    ZnLattice(shape; periodic = false)

Construct an `N`-dimensional cubic lattice of type `ZnLattice{N}`.

# Arguments
* `shape`: Tuple of sizes of all dimensions `(nx, ny, ...)`.
"""
function ZnLattice(shape::Tuple{Vararg{Int,N}};periodic::Bool = false) where N
    @assert N > 0 "The lattice must be positive-dimensional."
    if N == 1
        return ZnLattice{1}(Grid(collect(shape);periodic=periodic),shape,([1],),([shape[end]],),periodic)
    else
        L = Grid(collect(shape);periodic=periodic);
        lids = LinearIndices(shape)
        symbids = [1;fill(:(:), N-1)]
        Vint = Tuple([Core.eval(CornerSpaceRenorm,Expr(:ref, :($lids), circshift(symbids,d)...))[:] for d in 0:N-1])
        symbids = [:end;fill(:(:), N-1)]
        Vext = Tuple([Core.eval(Main,Expr(:ref, :($lids), circshift(symbids,d)...))[:] for d in 0:N-1])

        return ZnLattice{N}(L,shape,Vint,Vext,periodic)
    end
end

# Extend AbstractGraph{Int64} methods to all Lattice subtypes.
@inline Base.eltype(L::Lattice) = Base.eltype(L.L);
@inline LightGraphs.vertices(L::Lattice) = LightGraphs.vertices(L.L);
@inline LightGraphs.edges(L::Lattice) = LightGraphs.edges(L.L);
@inline LightGraphs.nv(L::Lattice) = LightGraphs.nv(L.L);
@inline LightGraphs.ne(L::Lattice) = LightGraphs.ne(L.L);
@inline LightGraphs.has_edge(L::Lattice, s::Int64, d::Int64) = LightGraphs.has_edge(L.L,s,d);
@inline LightGraphs.has_vertex(L::Lattice) = LightGraphs.has_vertex(L.L);
@inline LightGraphs.has_edge(L::Lattice, e::AbstractEdge{Int64}) = LightGraphs.has_edge(L.L,e);
@inline LightGraphs.inneighbors(L::Lattice) = LightGraphs.inneighbors(L.L);
@inline LightGraphs.has_vertex(L::Lattice, v::Int64) = LightGraphs.has_vertex(L.L,v);
@inline LightGraphs.outneighbors(L::Lattice) = LightGraphs.outneighbors(L.L);
@inline LightGraphs.inneighbors(L::Lattice, v::Int64) = LightGraphs.inneighbors(L.L,v);
@inline LightGraphs.outneighbors(L::Lattice, v::Int64) = LightGraphs.outneighbors(L.L,v);
@inline LightGraphs.is_directed(L::Lattice) = false;
@inline LightGraphs.is_directed(T::Type{L}) where {L<:Lattice} = false;

GraphPlot.gplot(L::Lattice; kwargs...) = GraphPlot.gplot(L.L; kwargs...)
GraphPlot.gplot(L::Lattice, locs_x_in::Vector{R}, locs_y_in::Vector{R}; kwargs...) where R<:Real = GraphPlot.gplot(L.L, locs_x_in, locs_y_in; kwargs...)

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

function rem_boundaries!(L::ZnLattice{N}, d::Int) where N
    if L.pbc == true && L.shape[d] > 2
        eb = [e for e in edges(L) if e.src ∈ L.Vint[d] && e.dst ∈ L.Vext[d]]
        for e in eb
            rem_edge!(L.L,e)
        end
    end
    return L
end

function add_boundaries!(L::ZnLattice{N}, d::Int) where N
    if L.pbc == true
        for i in 1:length(L.Vext[d])
            add_edge!(L.L,Edge(L.Vext[d][i],L.Vint[d][i]))
        end
    end
    return L
end
"""
    union(L1, L2, d)

Merge two `ZnLattice{N}`s along the `d`-th dimension.
"""
function Base.union(L1::ZnLattice{N},L2::ZnLattice{N},d::Int) where N
    @assert d <= N "Merging direction should be smaller than that of the lattices."
    @assert length(L1.Vext[d]) == length(L2.Vint[d]) "Lattices cannot be merged in this direction."
    @assert L1.pbc == L2.pbc "Cannot merge a periodic and an open lattices together."
    L1 = deepcopy(L1)
    L2 = deepcopy(L2)
    # Remove boundary edges along the merging dimension if periodic
    rem_boundaries!(L1, d)
    rem_boundaries!(L2, d)
    pbc = L1.pbc

    L = blockdiag(L1.L,L2.L)
    for i in 1:length(L1.Vext[d])
        add_edge!(L,Edge(L1.Vext[d][i],nv(L1)+L2.Vint[d][i]))
    end
    shape = collect(L1.shape);
    shape[d] += L2.shape[d]
    d⊥ = [_d for _d in 1:N if _d!= d]

    Vint = collect(L1.Vint)
    Vext = [nv(L1) .+ vs for vs in L2.Vext]
    for _d in d⊥
        Vint[_d] = [L1.Vint[_d]; nv(L1) .+ L2.Vint[_d]]
        Vext[_d] = [L1.Vext[_d]; nv(L1) .+ L2.Vext[_d]]
    end

    lat = ZnLattice{N}(L,Tuple(shape),Tuple(Vint),Tuple(Vext),pbc)
    # Add edges between boundaries if periodic
    add_boundaries!(lat, d)
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
              O3<:AbstractOperator{B,B},O4<:AbstractOperator{B,B}} <: AbstractSystem
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
    observables::Vector{Dict{String,O4}}
end

"""
    System(lat, H, lHt, J)

Construct a `System` from a `Lattice` and operators defining the model.
# Arguments
* `lat`: `Lattice`.
* `H`: Hamiltonian of the full system.
* `lHt`: local tunneling operator.
* `J`: jump operators of the full system.
"""
function System(lat::L,H::O1,lHt::O2,J::Vector{O3},obs::Union{Vector{Dict{String,O4}},Missing}=missing) where {N,L<:Lattice,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B},
                                                       O4<:AbstractOperator{B,B}}
    gbasis = H.basis_l;
    @assert nv(lat) == length(gbasis.bases)
    @assert lHt.basis_l == first(gbasis.bases)
    @assert H.basis_r == H.basis_l == gbasis
    @assert first(J).basis_l == gbasis
    @assert first(J).basis_r == gbasis
    @assert ismissing(obs) || length(obs) == nv(lat)
    Httop    = [embed(gbasis,[i],[lHt]) for i in lat.Vtop]
    Htbottom = [embed(gbasis,[i],[lHt]) for i in lat.Vbottom]
    Htleft   = [embed(gbasis,[i],[lHt]) for i in lat.Vleft]
    Htright  = [embed(gbasis,[i],[lHt]) for i in lat.Vright]
    if ismissing(obs)
        return System{L,B,O1,eltype(Httop),O3,typeof(H)}(lat,gbasis,H,Httop,Htbottom,Htleft,Htright,J,[Dict{String,typeof(H)}() for i in 1:nv(lat)])
    else
        return System{L,B,O1,eltype(Httop),O3,O4}(lat,gbasis,H,Httop,Htbottom,Htleft,Htright,J,obs)
    end
end

"""
    ZnSystem <: AbstractSystem

Quantum dissipative system defined on a `ZnLattice`.
"""
struct ZnSystem{N,L<:ZnLattice{N},B<:Basis,
                O1<:AbstractOperator{B,B},O2<:AbstractOperator{B,B},
                O3<:AbstractOperator{B,B},O4<:AbstractOperator{B,B}} <: AbstractSystem
    # Lattice
    lattice::L
    # Bases
    gbasis::B
    # Operators
    H::O1
    trate::Tuple{Vararg{ComplexF64,N}}
    Htint::Tuple{Vararg{Vector{O2},N}}
    Htext::Tuple{Vararg{Vector{O2},N}}
    J::Vector{O3}
    observables::Vector{Dict{String,O4}}
end

"""
    ZnSystem(lat, H, lHt, J, obs)

Construct a `ZnSystem` from a `ZnLattice` and operators defining the model.
# Arguments
* `lat`: `ZnLattice`.
* `H`: Hamiltonian of the full system.
* `lHt`: local tunneling operator.
* `J`: jump operators of the full system.
* `obs`: Observables to be transformed to the corner space together with the
Liouvillian. Can be a `Dict{String,<local operator type>}` of local operators or
an array of `nv(lat)` `Dict{String,<global operator type>}`s of operators in the
global basis.
"""
function ZnSystem(lat::ZnLattice{M},H::O1,trate::Tuple{Vararg{Number,M}},lHt::O2,J::Vector{O3},obs::Union{Vector{Dict{String,O4}},Missing}=missing) where {M,N,L<:Lattice,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B},
                                                       O4<:AbstractOperator{B,B}}
    gbasis = H.basis_l;
    @assert nv(lat) == length(gbasis.bases)
    @assert lHt.basis_l == first(gbasis.bases)
    @assert H.basis_r == H.basis_l == gbasis
    @assert first(J).basis_l == gbasis
    @assert first(J).basis_r == gbasis
    _trate::Tuple{Vararg{ComplexF64,M}} = Tuple(ComplexF64.(collect(trate)))
    Htint = Tuple([[embed(gbasis,v,lHt) for v in lat.Vint[d]] for d in 1:M])
    Htext = Tuple([[embed(gbasis,v,lHt) for v in lat.Vext[d]] for d in 1:M])

    if ismissing(obs)
        return ZnSystem{M,ZnLattice{M},B,O1,eltype(first(Htint)),O3,typeof(H)}(lat,gbasis,H,_trate,Htint,Htext,J,[Dict{String,typeof(H)}() for i in 1:nv(lat)])
    else
        return ZnSystem{M,ZnLattice{M},B,O1,eltype(first(Htint)),O3,O4}(lat,gbasis,H,_trate,Htint,Htext,J,obs)
    end
end

function ZnSystem(lat::ZnLattice{M},H::O1,trate::Tuple{Vararg{Number,M}},lHt::O2,J::Vector{O3},lobs::Dict{String,O4}) where {M,N,L<:Lattice,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B},
                                                       O4<:AbstractOperator{LB,LB}}
    gbasis = H.basis_l;
    obs = [Dict([name => embed(gbasis,i,lop) for (name,lop) in lobs]) for i in 1:nv(lat)]
    return ZnSystem(lat,H,Tuple(ComplexF64.(collect(trate))),lHt,J,obs)
end

GraphPlot.gplot(s::AbstractSystem; kwargs...) = GraphPlot.gplot(s.lattice.L; kwargs...)
GraphPlot.gplot(s::AbstractSystem, locs_x_in::Vector{R}, locs_y_in::Vector{R}; kwargs...) where R<:Real = GraphPlot.gplot(s.lattice.L, locs_x_in, locs_y_in; kwargs...)

"""
    plot_system(s; kwargs...)

Plot the lattice on which `s` is supported, equivalent to `GraphPlot.gplot(s.lattice.L; kwargs...)`.
# Arguments
See documentation of [`GraphPlot.gplot`](@ref).
"""
plot_system(s::AbstractSystem; kwargs...) = GraphPlot.gplot(s.lattice.L; kwargs...)

"""
    plot_system(s, locs_x_in, locs_y_in; kwargs...)

Equivalent to `GraphPlot.gplot(s.lattice.L, locs_x_in, locs_y_in; kwargs...)`.
# Arguments
See documentation of [`GraphPlot.gplot`](@ref).
"""
plot_system(s::AbstractSystem, locs_x_in::Vector{R}, locs_y_in::Vector{R}; kwargs...) where R<:Real = GraphPlot.gplot(s.lattice.L, locs_x_in, locs_y_in; kwargs...)
