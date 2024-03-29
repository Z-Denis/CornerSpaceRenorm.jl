"""
    hamiltonian(L, lH, lHt)

Construct a Hamiltonian from a `Lattice` and operators.
# Arguments
* `L`: Lattice.
* `lH`: local Hamiltonian.
* `lHt`: `Tuple` of the form (t,Op) from which tunnelling terms are built
according to Σ<i,j> (t * Op[i] ⊗ Op[j]' + t' * Op[j] ⊗ Op[i]')
"""
function hamiltonian(L::Lattice, lH::O1, lHt::Tuple{Number,O2}) where {B<:Basis,O1<:AbstractOperator{B,B},O2<:AbstractOperator{B,B}}
    gbasis = CompositeBasis([lH.basis_l for i in 1:nv(L)]...);
    H = begin
            if !isdense(lH) && !isdense(lHt[2])
                SparseOperator(gbasis);
            else
                DenseOperator(gbasis);
            end
        end

    for e in edges(L)
        H.data .+= lHt[1] * embed(gbasis,[e.src, e.dst],[lHt[2]; dagger(lHt[2])]).data;
    end
    H.data .= H.data + H.data';

    for v in vertices(L)
        H.data .+= embed(gbasis, v, lH).data;
    end

    return H;
end
hamiltonian(L::NdLattice{N}, lH::O1, lHt::Tuple{Number,O2}) where {N,B<:Basis,O1<:AbstractOperator{B,B},O2<:AbstractOperator{B,B}} = hamiltonian(L::NdLattice{N}, lH::O1, lHt[1], lHt[2])


#function hamiltonian(L::Lattice, lH::O1, trate::Number, lHt::O2) where {N,Q,B<:Basis,O1<:AbstractOperator{B,B},O2<:AbstractOperator{B,B}}
#    hamiltonian(L, lH, (trate, lHt))
#end

"""
    hamiltonian(L, lH, trate, lHt)

Construct a Hamiltonian from a `Lattice` and operators.
# Arguments
* `L`: `NdLattice{N}`.
* `lH`: local Hamiltonian.
* `trate`: tunnelling rate or N-`Tuple` containing a rate per lattice dimension.
* `lHt`: local operator from which tunnelling terms are built as lHt[i]' ⊗ lHt[j]'.
"""
function hamiltonian(L::NdLattice{N}, lH::O1, trate::Tuple{Vararg{Number,Q}}, lHt::O2) where {N,Q,B<:Basis,O1<:AbstractOperator{B,B},O2<:AbstractOperator{B,B}}
    @assert Q == N "Number of tunnelling rates must match number of dimensions"
    gbasis = CompositeBasis([lH.basis_l for i in 1:nv(L)]...);
    H = begin
            if !isdense(lH) && !isdense(lHt)
                SparseOperator(gbasis);
            else
                DenseOperator(gbasis);
            end
        end
    lids = LinearIndices(L.shape)
    function lattice_slice(d::Int,i::Int)
        symbids = [mod1(i,L.shape[d]);fill(:(:), N-1)]
        Core.eval(CornerSpaceRenorm,Expr(:ref, :($lids), circshift(symbids,d-1)...))[:]
    end
    # E[d]: set of edges along the d-th dimension
    E = [Vector{LightGraphs.SimpleGraphs.SimpleEdge{Int64}}(undef,0) for i in 1:N]
    for d in 1:N
        i_max = L.pbc && (L.shape[d] > 2) ? L.shape[d] : L.shape[d]-1
        E[d] = vcat([[Edge(e) for e in zip(lattice_slice(d,i),lattice_slice(d,i+1))] for i in 1:i_max]...)
        if L.pbc
            if L.shape[d] == 1
                loops = [Edge(v=>v) for v in vertices(L)]
                E[d] = [loops; loops]
            elseif L.shape[d] == 2
                E[d] = [E[d]; E[d]]
            end
        end
    end

    for d in 1:N
        for e in E[d]
            H.data .+= trate[d] * embed_safe(gbasis, e, [lHt; dagger(lHt)]).data;
        end
    end
    H.data .= H.data + H.data';

    for v in vertices(L)
        H.data .+= embed(gbasis, v, lH).data;
    end

    return H;
end

function hamiltonian(L::NdLattice{N}, lH::O1, trate, lHt) where {N,B<:Basis,O1<:AbstractOperator{B,B}}
    gbasis = CompositeBasis([lH.basis_l for i in 1:nv(L)]...);

    _trate = begin
        if typeof(trate) <: Number
            Tuple(fill(ComplexF64(trate), N))
        elseif typeof(trate)<:Tuple && eltype(trate)<:Number
            @assert length(trate) == N "Number of coupling rates must match number of dimensions"
            ComplexF64.(trate)
        else
            error("The coupling rate must either be a Number or a Tuple")
        end
    end

    _lHt = begin
        if typeof(lHt) <: AbstractOperator{B,B}
            Tuple([lHt] for d in 1:N)
        elseif typeof(lHt) <: Vector{O} where {O<:AbstractOperator{B,B}}
            Tuple(deepcopy(lHt) for d in 1:N)
        elseif typeof(lHt)<:Tuple
            @assert length(lHt) == N "Number of coupling operators must match number of dimensions"
            if eltype(lHt)<:AbstractOperator{B,B}
                Tuple([lHt[d]] for d in 1:N)
            elseif eltype(lHt)<:Vector{O} where {O<:AbstractOperator{B,B}}
                Tuple([h for h in lHt[d]] for d in 1:N)
            else
                error("The tuple of coupling operators must contain either local operators or a vector or local operators")
            end
        else
            error("The coupling operator must be a local operator, a tuple of local operators or a tuple of vectors of local operators")
        end
    end

    H = begin
            if !isdense(lH) && !isdense(first(first(_lHt)))
                SparseOperator(gbasis);
            else
                DenseOperator(gbasis);
            end
        end

    lids = LinearIndices(L.shape)
    function lattice_slice(d::Int,i::Int)
        symbids = [mod1(i,L.shape[d]);fill(:(:), N-1)]
        collect(Core.eval(CornerSpaceRenorm,Expr(:ref, :($lids), circshift(symbids,d-1)...)))[:]
    end
    # E[d]: set of edges along the d-th dimension
    E = [Vector{LightGraphs.SimpleGraphs.SimpleEdge{Int64}}(undef,0) for i in 1:N]
    for d in 1:N
        i_max = L.pbc && (L.shape[d] > 2) ? L.shape[d] : L.shape[d]-1
        E[d] = vcat([[Edge(e) for e in zip(lattice_slice(d,i),lattice_slice(d,i+1))] for i in 1:i_max]...)
        if L.pbc
            if L.shape[d] == 1
                loops = [Edge(v=>v) for v in vertices(L)]
                E[d] = [loops; loops]
            elseif L.shape[d] == 2
                E[d] = [E[d]; E[d]]
            end
        end
    end

    for d in 1:N,  e in E[d], h in _lHt[d]
        H.data .+= _trate[d] * embed_safe(gbasis, e, [h, dagger(h)]).data;
    end
    H.data .= H.data + H.data';

    for v in vertices(L)
        H.data .+= embed(gbasis, v, lH).data;
    end

    return H
end

function embed_safe(gb::CompositeBasis, e::AbstractEdge, ops::AbstractVector{O}) where O<:AbstractOperator
    if e.src != e.dst
        embed(gb,[e.src, e.dst],ops)
    else
        embed(gb,e.src,prod(ops))
    end
end

"""
    dissipators(L, J)

Construct the jump operators of a system from a `Lattice` and local jump operators.
# Arguments
* `L`: Lattice.
* `J`: Array of on-site local jump operators.
"""
function dissipators(L::Lattice, J::Vector{O}) where {B<:Basis,O<:AbstractOperator{B,B}}
    gbasis = CompositeBasis([first(J).basis_l for i in 1:nv(L)]...);
    return vcat([[embed(gbasis, v, J[i]) for i in 1:length(J)] for v in vertices(L)]...)
end

isdense(op::DataOperator) = typeof(op.data) <: DenseMatrix
isdense(op::AbstractOperator) = false
