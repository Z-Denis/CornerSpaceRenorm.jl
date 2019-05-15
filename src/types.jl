using LightGraphs, GraphPlot

const AbstractLattice = AbstractGraph{Int64}
"""
    Lattice <: AbstractLattice

Abstract supertype for all lattices.
"""
abstract type Lattice <: AbstractLattice end

"""
    NdLattice{N} <: Lattice

Z^`N` lattice type. Generic lattice with `N` directions.
"""
struct NdLattice{N} <: Lattice
    L::SimpleGraph{Int64}
    shape::Tuple{Vararg{Int,N}}

    Vint::Tuple{Vararg{Vector{Int64},N}}
    Vext::Tuple{Vararg{Vector{Int64},N}}

    # Boundary conditions
    pbc::Bool
end

"""
    NdLattice(shape; periodic = false)

Construct an `N`-dimensional cubic lattice of type `NdLattice{N}`.

# Arguments
* `shape`: Tuple of sizes of all dimensions `(nx, ny, ...)`.
"""
function NdLattice(shape::Tuple{Vararg{Int,N}};periodic::Bool = false) where N
    @assert N > 0 "The lattice must be positive-dimensional."
    if N == 1
        return NdLattice{1}(Grid(collect(shape);periodic=periodic),shape,([1],),([shape[end]],),periodic)
    else
        L = Grid(collect(shape);periodic=periodic);
        lids = LinearIndices(shape)
        symbids = [1;fill(:(:), N-1)]
        Vint = Tuple([Core.eval(CornerSpaceRenorm,Expr(:ref, :($lids), circshift(symbids,d)...))[:] for d in 0:N-1])
        symbids = [:end;fill(:(:), N-1)]
        Vext = Tuple([Core.eval(Main,Expr(:ref, :($lids), circshift(symbids,d)...))[:] for d in 0:N-1])

        return NdLattice{N}(L,shape,Vint,Vext,periodic)
    end
end

Base.show(io::IO, L::NdLattice) = L.pbc ? print(io,join([string(s) for s in L.shape], " x "), " NdLattice with PBC") : print(io,join([string(s) for s in L.shape], " x "), " NdLattice with OBC");
function Base.show(io::IO, ::MIME"text/plain", L::NdLattice)
    print(io,string(typeof(L)),":\n")
    print(io,"  Shape: ",join([string(s) for s in L.shape], " x "),"\n")
    print(io,"  Boundaries:\n")
    for d in 1:length(L.shape)
        print(io,"  d = ",d,":\n")
        print(io,"    Input vertices: ",L.Vint[d],"\n")
        print(io,"    Input vertices: ",L.Vext[d],"\n")
    end

    if L.pbc
        print(io,"  Boundary conditions: periodic")
    else
        print(io,"  Boundary conditions: open")
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::NdLattice) = x

# Extend AbstractGraph{Int64} methods to all Lattice subtypes.
@inline Base.eltype(L::Lattice) = Base.eltype(L.L);
@inline LightGraphs.vertices(L::Lattice) = LightGraphs.vertices(L.L);
@inline LightGraphs.edges(L::Lattice) = LightGraphs.edges(L.L);
@inline LightGraphs.nv(L::Lattice) = LightGraphs.nv(L.L);
@inline LightGraphs.ne(L::Lattice) = LightGraphs.ne(L.L);
@inline LightGraphs.has_edge(L::Lattice, s::Int64, d::Int64) = LightGraphs.has_edge(L.L,s,d);
@inline LightGraphs.has_edge(L::Lattice, e::AbstractEdge{Int64}) = LightGraphs.has_edge(L.L,e);
@inline LightGraphs.has_vertex(L::Lattice, v::Int64) = LightGraphs.has_vertex(L.L,v);
@inline LightGraphs.inneighbors(L::Lattice, v::Int64) = LightGraphs.inneighbors(L.L,v);
@inline LightGraphs.outneighbors(L::Lattice, v::Int64) = LightGraphs.outneighbors(L.L,v);
@inline LightGraphs.is_directed(L::Lattice) = false;
@inline LightGraphs.is_directed(T::Type{L}) where {L<:Lattice} = false;

GraphPlot.gplot(L::Lattice; kwargs...) = GraphPlot.gplot(L.L; kwargs...)
GraphPlot.gplot(L::Lattice, locs_x_in::Vector{R}, locs_y_in::Vector{R}; kwargs...) where R<:Real = GraphPlot.gplot(L.L, locs_x_in, locs_y_in; kwargs...)

function rem_boundaries!(L::NdLattice{N}, d::Int) where N
    if L.pbc == true && L.shape[d] > 2
        eb = [e for e in edges(L) if e.src ∈ L.Vint[d] && e.dst ∈ L.Vext[d]]
        for e in eb
            rem_edge!(L.L,e)
        end
    end
    return L
end

function add_boundaries!(L::NdLattice{N}, d::Int) where N
    if L.pbc == true
        for i in 1:length(L.Vext[d])
            add_edge!(L.L,Edge(L.Vext[d][i],L.Vint[d][i]))
        end
    end
    return L
end
"""
    union(L1, L2, d)

Merge two `NdLattice{N}`s along the `d`-th dimension.
"""
function Base.union(L1::NdLattice{N},L2::NdLattice{N},d::Int) where N
    @assert d <= N "Merging direction should be smaller than that of the lattices."
    length(L1.Vext[d]) == length(L2.Vint[d]) || throw(DimensionMismatch("Lattices cannot be merged in this direction."))
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

    lat = NdLattice{N}(L,Tuple(shape),Tuple(Vint),Tuple(Vext),pbc)
    # Add edges between boundaries if periodic
    add_boundaries!(lat, d)
end

"""
    AbstractSystem

Abstract supertype for all systems.
"""
abstract type AbstractSystem end

"""
    NdSystem <: AbstractSystem

Quantum dissipative system defined on a `NdLattice`.
"""
struct NdSystem{N,L<:NdLattice{N},B<:Basis,
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
    NdSystem(lat, H, trate, lHt, J, obs)

Construct a `NdSystem` from a `NdLattice` and operators defining the model.
# Arguments
* `lat`: `NdLattice{N}`.
* `H`: Hamiltonian of the full system.
* `trate`: tunneling rate or N-`Tuple` of tunneling rates along each dimension.
* `lHt`: local tunneling operator.
* `J`: jump operators of the full system.
* `obs`: Observables to be transformed to the corner space together with the
Liouvillian. Can be a `Dict{String,<local operator type>}` of local operators or
an array of `nv(lat)` `Dict{String,<global operator type>}`s of operators in the
global basis.
"""
function NdSystem(lat::NdLattice{M},H::O1,trate::Tuple{Vararg{Number,Q}},lHt::O2,J::Vector{O3},obs::Union{Vector{Dict{String,O4}},Missing}=missing) where {M,N,Q,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B},
                                                       O4<:AbstractOperator{B,B}}
    gbasis = H.basis_l;
    @assert nv(lat) == length(gbasis.bases)
    @assert Q == M "Number of tunnelling rates must match number of dimensions"
    @assert lHt.basis_l == first(gbasis.bases) "Global basis of H must match the local basis of lHt"
    @assert H.basis_r == H.basis_l == gbasis
    @assert first(J).basis_l == first(J).basis_r == gbasis
    _trate::Tuple{Vararg{ComplexF64,M}} = Tuple(ComplexF64.(collect(trate)))
    Htint = Tuple([[embed(gbasis,v,lHt) for v in lat.Vint[d]] for d in 1:M])
    Htext = Tuple([[embed(gbasis,v,lHt) for v in lat.Vext[d]] for d in 1:M])

    if ismissing(obs)
        return NdSystem{M,NdLattice{M},B,O1,eltype(first(Htint)),O3,typeof(H)}(lat,gbasis,H,_trate,Htint,Htext,J,[Dict{String,typeof(H)}() for i in 1:nv(lat)])
    else
        return NdSystem{M,NdLattice{M},B,O1,eltype(first(Htint)),O3,O4}(lat,gbasis,H,_trate,Htint,Htext,J,obs)
    end
end

function NdSystem(lat::NdLattice{M},H::O1,trate::Number,lHt::O2,J::Vector{O3},obs::Union{Vector{Dict{String,O4}},Missing}=missing) where {M,N,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B},
                                                       O4<:AbstractOperator{B,B}}
    return NdSystem(lat,H,Tuple([trate for i in 1:M]),lHt,J,obs)
end

function NdSystem(lat::NdLattice{M},H::O1,trate::Tuple{Vararg{Number,Q}},lHt::O2,J::Vector{O3},lobs::Dict{String,O4}) where {M,N,Q,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B},
                                                       O4<:AbstractOperator{LB,LB}}
    gbasis = H.basis_l;
    @assert Q == M "Number of rates must match number of dimensions."
    obs = [Dict([name => embed(gbasis,i,lop) for (name,lop) in lobs]) for i in 1:nv(lat)]
    return NdSystem(lat,H,Tuple(ComplexF64.(collect(trate))),lHt,J,obs)
end

function NdSystem(lat::NdLattice{M},H::O1,trate::Number,lHt::O2,J::Vector{O3},lobs::Dict{String,O4}) where {M,N,LB<:Basis,
                                                       B<:CompositeBasis{Tuple{Vararg{LB,N}}},
                                                       O1<:AbstractOperator{B,B},
                                                       O2<:AbstractOperator{LB,LB},
                                                       O3<:AbstractOperator{B,B},
                                                       O4<:AbstractOperator{LB,LB}}
    return NdSystem(lat,H,Tuple([trate for i in 1:M]),lHt,J,lobs)
end

function Base.show(io::IO, s::NdSystem)
    s.lattice.pbc ? print(io,join([string(_s) for _s in s.lattice.shape], " x "), " NdSystem with PBC") : print(io,join([string(_s) for _s in s.lattice.shape], " x "), " NdSystem with OBC");
    if typeof(s.gbasis) <: CompositeBasis && all([s.gbasis.bases[i] == s.gbasis.bases[1] for i in 1:length(s.gbasis.bases)])
        print(io, " and local basis: ")
        Base.show(io, first(s.gbasis.bases))
    else
        print(io, " and global basis: ")
        Base.show(io, s.gbasis)
    end
end
function Base.show(io::IO, ::MIME"text/plain", s::NdSystem)
    print(io, "NdSystem:\n")
    print(io,"  Lattice: ")
    Base.show(io,s.lattice)
    if typeof(s.gbasis) <: CompositeBasis && all([s.gbasis.bases[i] == s.gbasis.bases[1] for i in 1:length(s.gbasis.bases)])
        print(io, "\n  Homogeneous local basis: ")
        Base.show(io, first(s.gbasis.bases))
    else
        print(io, "\n  Global basis: ")
        Base.show(io, s.gbasis)
    end
    print(io,"\n  Hamiltonian: ",join([string(_s) for _s in size(s.H)], " x ")," ",typeof(s.H).name,"\n")
    print(io,"  ",length(s.J)," lindblad ops: ",join([string(_s) for _s in size(first(s.J))], " x ")," ",typeof(first(s.J)).name,"\n")
    print(io, "  Tunneling rates: \n")
    for d in 1:length(s.lattice.shape)
        print(io,"    d = ",d,": ",s.trate[d],"\n")
    end
    obs_names = unique(keys.(s.observables))
    if length(obs_names) == 1
        print(io, "  Local observables: ",join([obs for obs in obs_names[1]], ", "))
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::NdSystem) = x

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

"""
    CornerBasis <: Basis

Corner basis.
"""
struct CornerBasis <: Basis
    shape::Vector{Int64}
    M::Int64
    local_shapes::Vector{Vector{Int64}}
end

"""
    CornerBasis(M)

Construct a corner basis with corner dimension `M`.
"""
CornerBasis(M::Int64) = CornerBasis([M], M, [[M]])
"""
    CornerBasis(b1, b2, M)

Construct a corner basis with corner dimension `M` from two basis.
"""
CornerBasis(b1::Basis, b2::Basis, M::Int64) = CornerBasis([M], M, [b1.shape, b2.shape])
CornerBasis(b1::CornerBasis, b2::CornerBasis, M::Int64) = CornerBasis([M], M, vcat(b1.local_shapes, b2.local_shapes))

import Base: ==
==(b1::CornerBasis, b2::CornerBasis) = (b1.M .== b2.M) &&
                                       (b1.shape == b2.shape) &&
                                       all(b1.local_shapes .== b2.local_shapes)

Base.show(io::IO, cb::CornerBasis) = print(io,"CornerBasis with dimension M = ",cb.M);
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::CornerBasis) = x
