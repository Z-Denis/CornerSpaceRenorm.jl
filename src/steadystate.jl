using IterativeSolvers, LinearMaps, LinearAlgebra

"""
    steadystate_bicg(H, J, l=2; [log=false], kwargs...) -> rho[, log]

Compute the steady state density matrix of a Hamiltonian and a set of jump
operators by solving `L rho = 0` via the stabilized biconjugate gradient method.
The first line of the Liouvillian is overwritten to enforce a non-trivial trace
one solution, this approximation yields an error of the order of the inverse of
the square of the size of the Hilbert space.

# Arguments
* `H`: dense Hamiltonian.
* `J`: array of dense jump operators.
* `l`: number of GMRES steps per iteration.
* `kwargs...`: Further arguments are passed on to the iterative solver.

See also: [`CornerSpaceRenorm.steadystate_bicg_LtL`](@ref)
"""
function steadystate_bicg(H::AbstractOperator{B,B}, J::Vector{O}, l::Int=2; log::Bool=false, kwargs...) where {B<:Basis,O<:AbstractOperator{B,B}}
    @assert isdense(H) "The Hamiltonian must be dense"
    @assert all(isdense.(J)) "The jump operators must be dense"
    ρ0 = DenseOperator(H.basis_l)
    ρ0.data[1,1] = ComplexF64(1.0)
    return steadystate_bicg!(ρ0,H,J,l;log=log,kwargs...)
end

"""
    steadystate_bicg(s, l=2; [log=false], kwargs...) -> rho[, log]

Compute the steady state density matrix of an `AbstractSystem` by solving
`L rho = 0` via the stabilized biconjugate gradient method.
The first line of the Liouvillian is overwritten to enforce a non-trivial trace
one solution, this approximation yields an error of the order of the inverse of
the square of the size of the Hilbert space.

# Arguments
* `s`: Instance of any subtype of `AbstractSystem`.
* `l`: number of GMRES steps per iteration.
* `kwargs...`: Further arguments are passed on to the iterative solver.

See also: [`CornerSpaceRenorm.steadystate_bicg_LtL`](@ref)
"""
steadystate_bicg(s::AbstractSystem, l::Int=2; log::Bool=false, kwargs...) = steadystate_bicg(DenseOperator(s.H), DenseOperator.(s.J), l; log=log, kwargs...)

"""
    steadystate_bicg!(rho0, H, J, l=2; [log=false], kwargs...) -> rho[, log]

Compute the steady state density matrix of a Hamiltonian and a set of jump
operators by solving `L rho = 0` via the stabilized biconjugate gradient method.
The first line of the Liouvillian is overwritten to enforce a non-trivial trace
one solution, this approximation yields an error of the order of the inverse of
the square of the size of the Hilbert space.

# Arguments
* `rho0`: Initial guess.
* `H`: dense Hamiltonian.
* `J`: array of dense jump operators.
* `l`: number of GMRES steps per iteration.
* `kwargs...`: Further arguments are passed on to the iterative solver.

See also: [`steadystate_bicg`](@ref)
"""
function steadystate_bicg!(ρ0::AbstractOperator{B,B}, H::AbstractOperator{B,B}, J::Vector{O}, l::Int=2; log::Bool=false, tol::Float64 = sqrt(eps(real(ComplexF64))), kwargs...) where {B<:Basis,O<:AbstractOperator{B,B}}
    @assert isdense(ρ0) "The initial density matrix must be dense"
    @assert isdense(H) "The Hamiltonian must be dense"
    @assert all(isdense.(J)) "The jump operators must be dense"
    # Size of the Hilbert space
    M::Int = size(H.data,1)
    # Non-Hermitian Hamiltonian
    iHnh::Matrix{ComplexF64} = -im*H.data
    for i in 1:length(J)
        iHnh .+= -0.5adjoint(J[i].data)*J[i].data
    end

    # In-place update of y = Lx where L and x are respectively the vectorized
    # Liouvillian and the vectorized density matrix. y[1] is set to the trace
    # of the density matrix so as to enforce a trace one non-trivial solution.
    function mvecmul!(y::AbstractVector, x::AbstractVector)
        y .= zero(eltype(x));
        ym::Matrix{ComplexF64} = reshape(y,M,M)
        ρ::Matrix{ComplexF64} = reshape(x,M,M)
        Jρ_cache::Matrix{ComplexF64} = similar(ρ)

        BLAS.gemm!('N', 'N', one(ComplexF64), iHnh, ρ, one(ComplexF64), ym)
        BLAS.gemm!('N', 'C', one(ComplexF64), ρ, iHnh, one(ComplexF64), ym)
        @inbounds @views for i in 1:length(J)
            BLAS.gemm!('N','N', one(ComplexF64), J[i].data, ρ, zero(ComplexF64), Jρ_cache)
            BLAS.gemm!('N','C', one(ComplexF64), Jρ_cache, J[i].data, one(ComplexF64), ym)
        end
        y .= reshape(ym,size(y))
        y[1] = tr(ρ)

        return y
    end
    # Solution x must satisfy L.x = y with y[1] = tr(x) = 1 and y[j≠1] = 0.
    x0::Vector{ComplexF64} = reshape(ρ0.data,M^2)
    y::Vector{ComplexF64}  = zeros(ComplexF64,M^2)
    y[1] = one(ComplexF64)

    # Define the linear map lm: ρ ↦ L(ρ)
    lm = LinearMap{eltype(H.data)}(mvecmul!, length(y)::Int, length(y)::Int; ismutating=true, issymmetric=false, ishermitian=false, isposdef=false)

    # Perform the stabilized biconjugate gradient procedure and devectorize ρ
    res0_norm::Float64 = norm(mvecmul!(similar(y),x0) .- y)
    tol /= res0_norm + eps(Float64)
    if !log
        ρ0.data .= reshape(bicgstabl!(x0,lm,y,l;tol=tol,kwargs...),(M,M))
        return ρ0
    else
        R::Vector{ComplexF64}, history = bicgstabl!(x0,lm,y,l;log=log,tol=tol,kwargs...)
        ρ0.data .= reshape(R,(M,M))
        return ρ0, history
    end
end

"""
    steadystate_bicg!(rho0, s, l=2; [log=false], kwargs...) -> rho[, log]

Compute the steady state density matrix of an `AbstractSystem` by solving
`L rho = 0` via the stabilized biconjugate gradient method.
The first line of the Liouvillian is overwritten to enforce a non-trivial trace
one solution, this approximation yields an error of the order of the inverse of
the square of the size of the Hilbert space.

# Arguments
* `rho0`: Initial guess.
* `s`: Instance of any subtype of `AbstractSystem`.
* `l`: number of GMRES steps per iteration.
* `kwargs...`: Further arguments are passed on to the iterative solver.

See also: [`steadystate_bicg`](@ref)
"""
steadystate_bicg!(ρ0::AbstractOperator{B,B}, s::AbstractSystem, l::Int=2; log::Bool=false, tol::Float64 = sqrt(eps(real(ComplexF64))), kwargs...) where {B<:Basis} = steadystate_bicg!(ρ0, DenseOperator(s.H), DenseOperator.(s.J), l; log=log, tol=tol, kwargs...)

"""
    CornerSpaceRenorm.steadystate_bicg_LtL(H, J, l=2; [log=false], kwargs...) -> rho[, log]

Compute the steady state density matrix of a Hamiltonian and a set of jump
operators by solving `<L,L(rho)> + <tr,rho> tr = tr` via the stabilized
biconjugate gradient method. No approximation of the Liovillian is made in order
to enfore the trace one but convergence is slower and poorer.

# Arguments
* `H`: dense Hamiltonian.
* `J`: array of dense jump operators.
* `l`: number of GMRES steps per iteration.
* `kwargs...`: Further arguments are passed on to the iterative solver.

See also: [`steadystate_bicg`](@ref)
"""
function steadystate_bicg_LtL(H::AbstractOperator{B,B}, J::Vector{O}, l::Int=2; log::Bool=false, kwargs...) where {B<:Basis,O<:AbstractOperator{B,B}}
    @assert isdense(H) "The Hamiltonian must be dense"
    @assert all(isdense.(J)) "The jump operators must be dense"
    # Size of the Hilbert space
    M::Int = size(H.data,1)
    # Non-Hermitian Hamiltonian
    iHnh::Matrix{ComplexF64} = -im*H.data
    for i in 1:length(J)
        iHnh .+= -0.5adjoint(J[i].data)*J[i].data
    end

    T::Vector{ComplexF64} = zeros(ComplexF64,M^2)
    for i in 0:M-1
        T[(M+1)*i+1] = one(ComplexF64)
    end

    # In-place update of y = L.x where L and x are respectively the vectorized
    # Liouvillian and the vectorized density matrix. y[1] is set to the trace
    # of the density matrix so as to enforce a trace one non-trivial solution.
    function mvecmul!(y::AbstractVector, x::AbstractVector)
        y .= zero(eltype(x));
        ym::Matrix{ComplexF64} = reshape(y,M,M)
        ρ::Matrix{ComplexF64} = reshape(x,M,M)
        Jρ_cache::Matrix{ComplexF64} = similar(ρ)
        y .= tr(ρ) .* T

        BLAS.gemm!('N', 'N', one(ComplexF64), iHnh, ρ, one(ComplexF64), ym)
        BLAS.gemm!('N', 'C', one(ComplexF64), ρ, iHnh, one(ComplexF64), ym)
        @inbounds @views for i in 1:length(J)
            BLAS.gemm!('N','N', one(ComplexF64), J[i].data, ρ, zero(ComplexF64), Jρ_cache)
            BLAS.gemm!('N','C', one(ComplexF64), Jρ_cache, J[i].data, one(ComplexF64), ym)
        end

        ρ .= reshape(ym, size(ρ))

        BLAS.gemm!('C', 'N', one(ComplexF64), iHnh, ρ, one(ComplexF64), ym)
        BLAS.gemm!('N', 'N', one(ComplexF64), ρ, iHnh, one(ComplexF64), ym)
        @inbounds @views for i in 1:length(J)
            BLAS.gemm!('N','N', one(ComplexF64), ρ, J[i].data, zero(ComplexF64), Jρ_cache)
            BLAS.gemm!('C','N', one(ComplexF64), J[i].data, Jρ_cache, one(ComplexF64), ym)
        end

        y .+= reshape(ym,size(y))

        return y
    end
    # Solution x must satisfy L.x = y with y[1] = tr(x) = 1 and y[j≠1] = 0.
    y::Vector{ComplexF64} = zeros(ComplexF64,M^2)
    y .= T

    # Define the linear map lm: ρ ↦ L(ρ)
    lm = LinearMap{eltype(H.data)}(mvecmul!, length(y)::Int, length(y)::Int; ismutating=true, issymmetric=false, ishermitian=true, isposdef=true)

    # Perform the stabilized biconjugate gradient procedure and devectorize ρ
    if !log
        return DenseOperator(H.basis_l, reshape(bicgstabl(lm,y,l;kwargs...),(M,M)))
    else
        R::Vector{ComplexF64}, history = bicgstabl(lm,y,l;log=log,kwargs...)
        return DenseOperator(H.basis_l, reshape(R,(M,M))), history
    end
end


"""
    CornerSpaceRenorm.steadystate_bicg_LtL(s, l=2; [log=false], kwargs...) -> rho[, log]

Compute the steady state density matrix of an `AbstractSystem` by solving
`<L,L(rho)> + <tr,rho> tr = tr` via the stabilized biconjugate gradient method.
No approximation of the Liovillian is made in order to enfore the trace one but
convergence is slower and poorer.

# Arguments
* `s`: Instance of any subtype of `AbstractSystem`.
* `l`: number of GMRES steps per iteration.
* `kwargs...`: Further arguments are passed on to the iterative solver.

See also: [`steadystate_bicg`](@ref)
"""
steadystate_bicg_LtL(s::AbstractSystem, l::Int=2; log::Bool=false, kwargs...) = steadystate_bicg_LtL(DenseOperator(s.H), DenseOperator.(s.J), l; log=log, kwargs...)

"""
    steadystate.master(s; <keyword arguments>)

Calculate steady state using long time master equation evolution.

# Arguments
* `s`: Instance of any subtype of `AbstractSystem`.
* `kwargs...`: Further arguments are passed on to the iterative solver.
"""
QuantumOptics.steadystate.master(s::AbstractSystem; kwargs...) = QuantumOptics.steadystate.master(s.H, s.J; kwargs...)

#= Broken
function bloch_redfield_bicg(H::DenseOperator{B,B}, a_ops::Vector{O}, S::Function, l::Int=2; log::Bool=false, kwargs...) where {B<:Basis,O<:DenseOperator{B,B}}
    M::Int = size(H.data,1)

    e = eigen(Matrix(H.data))
    ω::Vector{Float64} = real.(e.values)
    V::Matrix{ComplexF64} = CornerSpaceRenorm.mgs(e.vectors)
    Sc::Matrix{Matrix{ComplexF64}} = [[S(i, j,_ω′ - _ω) for _ω in ω, _ω′ in ω] for i in 1:length(a_ops), j in 1:length(a_ops)]
    aH_ms::Vector{Matrix{ComplexF64}} = [V' * op.data * V for op in a_ops]
    bH_ms::Vector{Matrix{ComplexF64}} = [sum([Sc[j, i] .* aH_ms[j] for j in 1:length(a_ops)]) for i in 1:length(a_ops)]
    cH_ms::Vector{Matrix{ComplexF64}} = [sum([Sc[i, j] .* aH_ms[j] for j in 1:length(a_ops)]) for i in 1:length(a_ops)]

    function mvecmul!(y::AbstractVector, x::AbstractVector)
        y .= zero(eltype(x));
        ym::Matrix{ComplexF64} = reshape(y,M,M)
        ρH_m::Matrix{ComplexF64} = reshape(x,M,M)
        cache::Matrix{ComplexF64} = similar(ρH_m)

        @inbounds @views for i in 1:length(a_ops)
            BLAS.gemm!('N','N', one(ComplexF64), bH_ms[i], ρH_m, zero(ComplexF64), cache)
            BLAS.gemm!('N','N', one(ComplexF64), cache, aH_ms[i], one(ComplexF64), ym)

            BLAS.gemm!('N','N', -one(ComplexF64), aH_ms[i], cH_ms[i], zero(ComplexF64), cache)
            BLAS.gemm!('N','N',  one(ComplexF64), cache, ρH_m, one(ComplexF64), ym)
        end
        #BLAS.gemm!('N', 'N', -1.0im, H.data, ρH_m, one(ComplexF64), ym)
        ym .+= adjoint(ym)

        y .= reshape(ym,size(y))
        y[1] = tr(ρH_m)
        return y
    end
    # Solution x must satisfy L.x = y with y[1] = tr(x) = 1 and y[j≠1] = 0.
    y::Vector{ComplexF64} = zeros(ComplexF64,M^2)
    y[1] = one(ComplexF64)

    # Define the linear map lm: ρ ↦ L(ρ)
    lm = LinearMap{eltype(H.data)}(mvecmul!, length(y)::Int, length(y)::Int; ismutating=true, issymmetric=false, ishermitian=false, isposdef=false)

    # Perform the stabilized biconjugate gradient procedure and devectorize ρ
    if !log
        return DenseOperator(H.basis_l, V * reshape(bicgstabl(lm,y,l;kwargs...),(M,M)) * V')
    else
        R::Vector{ComplexF64}, history = bicgstabl(lm,y,l;log=log,kwargs...)
        return DenseOperator(H.basis_l, V * reshape(R,(M,M)) * V'), history
    end
end
=#
