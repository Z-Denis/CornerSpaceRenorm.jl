using CornerSpaceRenorm, QuantumOptics
using LightGraphs, GraphPlot
using LinearAlgebra
using Test, InteractiveUtils

@testset "CornerSpaceRenorm.jl" begin
    lb = FockBasis(1)
    L  = SquareLattice(2,2)
    gb = CompositeBasis([lb for i in 1:nv(L)])
    tH = (destroy(lb),create(lb))
    H  = hamiltonian(L,number(lb) + 1e-1 * (create(lb) + destroy(lb)),[tH,dagger.(tH)])
    J  = dissipators(L,[destroy(lb),1e-2 * create(lb)])

    # Test method extensions
    # TO DO: extend tests to instances of all subtypes of Lattice
    @test eltype(L) == eltype(L.L)
    @test vertices(L) == vertices(L.L)
    @test edges(L) == edges(L.L)
    @test nv(L) == nv(L.L)
    @test ne(L) == ne(L.L)
    @test all([has_edge(L,v1,v2) == has_edge(L.L,v1,v2) for v1 in vertices(L), v2 in vertices(L)])
    @test all([has_edge(L,e) == has_edge(L.L,e) for e in edges(L)])
    @test all([has_vertex(L,v) == has_vertex(L.L,v) for v in collect(vertices(L)) ∪ collect(-10:-1)])
    @test all([inneighbors(L,v) == inneighbors(L.L,v) for v in vertices(L)])
    @test all([outneighbors(L,v) == outneighbors(L.L,v) for v in vertices(L)])
    @test is_directed(L) == false
    @test all(is_directed.(subtypes(Lattice)) .== false)

    # TO DO: find a way to compare outputs
    #@test gplot(L) == gplot(L.L)
    #s = System(L,H,destroy(lb),J)
    #@test plot_system(s) == gplot(s.lattice.L)

    # Test merging
    function compare_systems(s1,s2)
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
        nothing
    end
    L1  = SquareLattice(2,1)
    L2  = SquareLattice(3,1)
    s1 = L1 ∪ L2
    s2 = vunion(L1,L2)
    compare_systems(s1,s2)
    L2 = SquareLattice(2,3)
    s1 = L1 ∪ L2
    s2 =  hunion(L1,L2)
    compare_systems(s1,s2)

    s1 = System(L,H,destroy(lb),J)
    s1 = vmerge(s1,s1)
    s1 = hmerge(s1,s1)

    s2 = System(L,H,destroy(lb),J)
    s2 = vmerge(s2,s2)
    s2 = hmerge(s2,s2)

    compare_systems(s1,s2)

    # Test cornerized merging
    # Dumb version
    s1 = System(L,H,destroy(lb),J)
    ρ1 = steadystate.master(s1.H,s1.J)[2][end]
    s1 = vmerge(s1,s1)
    cspace = corner_subspace(ρ1,ρ1,10)[1]
    s1 = cornerize(s1,cspace)
    ρ1 = steadystate.master(s1.H,s1.J)[2][end]
    s1 = hmerge(s1,s1)
    cspace = corner_subspace(ρ1,ρ1,10)[1]
    s1 = cornerize(s1,cspace)

    # Wise version
    s2 = System(L,H,destroy(lb),J)
    ρ2 = steadystate.master(s2.H,s2.J)[2][end]
    s2 = vmerge(s2,s2,ρ2,ρ2,10)
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

    # Test ZnLattice and ZnSystem
    # Check equivalence with SquareLattice
    ZnL = ZnLattice((2,2))
    s3 = ZnSystem(ZnL,H,destroy(lb),J)
    ρ3 = steadystate.master(s3.H,s3.J)[2][end]
    s3 = merge(s3,s3,1,ρ3,ρ3,10)
    ρ3 = steadystate.master(s3.H,s3.J)[2][end]
    s3 = merge(s3,s3,2,ρ3,ρ3,10)

    @test maximum(abs2.(s3.H.data .- s2.H.data)) ≈ 0
    @test maximum([maximum(abs2.(s3.J[i].data .- s2.J[i].data)) for i in 1:length(s3.J)]) ≈ 0
    @test maximum([maximum(abs2.(s3.Htint[1][i].data .- s2.Httop[i].data)) for i in 1:length(s3.Htint[1])]) ≈ 0
    @test maximum([maximum(abs2.(s3.Htint[2][i].data .- s2.Htleft[i].data)) for i in 1:length(s3.Htint[1])]) ≈ 0
    @test maximum([maximum(abs2.(s3.Htext[1][i].data .- s2.Htbottom[i].data)) for i in 1:length(s3.Htext[1])]) ≈ 0
    @test maximum([maximum(abs2.(s3.Htext[2][i].data .- s2.Htright[i].data)) for i in 1:length(s3.Htext[1])]) ≈ 0

    # Check steady states
    ρ3 = steadystate.master(s3.H,s3.J)[2][end]
    ρ2 = steadystate.master(s2.H,s2.J)[2][end]
    @test fidelity(ρ3,ρ2) ≈ 1.

    # Test hermitianize!
    O = randoperator(GenericBasis(300))
    O_copy = deepcopy(O)
    CornerSpaceRenorm.hermitianize!(O)
    @test ishermitian(O.data)
    @test maximum(abs2.(O.data .- (O_copy + dagger(O_copy)).data/2.)) ≈ 0

    # Test constructor with local observables
    O = randoperator(lb)
    lobs = Dict("O"=>O)
    s1 = ZnSystem(ZnL,H,destroy(lb),J,lobs)
    gobs = [Dict("O"=>embed(s1.gbasis, i, O)) for i in vertices(ZnL)]
    s2 = ZnSystem(ZnL,H,destroy(lb),J,gobs)
    @test s1.observables == s2.observables
end
