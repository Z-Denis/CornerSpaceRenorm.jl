using CornerSpaceRenorm, QuantumOptics
using Test

@testset "CornerSpaceRenorm.jl" begin
    lb = FockBasis(1)
    L  = SquareLattice(2,2)
    gb = CompositeBasis([lb for i in 1:nv(L)])
    tH = (destroy(lb),create(lb))
    H  = hamiltonian(L,number(lb) + 1e-1 * (create(lb) + destroy(lb)),[tH,dagger.(tH)])
    J  = dissipators(L,[destroy(lb),1e-2 * create(lb)])

    # Test merging
    s1 = System(L,H,destroy(lb),J)
    s1 = vmerge(s1,s1)
    s1 = hmerge(s1,s1)

    s2 = System(L,H,destroy(lb),J)
    s2 = vmerge(s2,s2)
    s2 = hmerge(s2,s2)

    @test typeof(s1) == typeof(s2)
    for f in fieldnames(typeof(s1))
        if f == :lattice
            L1 = getfield(s1,f)
            L2 = getfield(s2,f)
            for ff in fieldnames(typeof(L1))
                @test getfield(L1,ff) == getfield(L2,ff)
            end
        else
            @test getfield(s1,f) == getfield(s2,f)
        end
    end

    # Test cornerized merging
    # Dumb version
    s1 = System(L,H,destroy(lb),J)
    ρ1 = steadystate.master(s1.H,s1.J)[2][end]
    s1 = hmerge(s1,s1)
    cspace = corner_subspace(ρ1,ρ1,10)[1]
    s1 = cornerize(s1,cspace)
    ρ1 = steadystate.master(s1.H,s1.J)[2][end]
    s1 = hmerge(s1,s1)
    cspace = corner_subspace(ρ1,ρ1,10)[1]
    s1 = cornerize(s1,cspace)

    # Wise version
    s2 = System(L,H,destroy(lb),J)
    ρ2 = steadystate.master(s2.H,s2.J)[2][end]
    s2 = hmerge(s2,s2,ρ2,ρ2,10)
    ρ2 = steadystate.master(s2.H,s2.J)[2][end]
    s2 = hmerge(s2,s2,ρ2,ρ2,10)

    # Check that generated matrices are equivalent
    tol = 1e-20
    @test maximum(abs2.(s1.H.data .- s2.H.data)) < tol
    @test maximum([maximum(abs2.(s1.J[i].data .- s2.J[i].data)) for i in 1:length(s1.J)]) < tol
    @test maximum([maximum(abs2.(s1.Httop[i].data .- s2.Httop[i].data)) for i in 1:length(s1.Httop)]) < tol
    @test maximum([maximum(abs2.(s1.Htbottom[i].data .- s2.Htbottom[i].data)) for i in 1:length(s1.Htbottom)]) < tol
    @test maximum([maximum(abs2.(s1.Htleft[i].data .- s2.Htleft[i].data)) for i in 1:length(s1.Htleft)]) < tol
    @test maximum([maximum(abs2.(s1.Htright[i].data .- s2.Htright[i].data)) for i in 1:length(s1.Htright)]) < tol
end
