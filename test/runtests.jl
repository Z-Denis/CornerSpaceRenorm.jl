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
    s1 = CornerSpaceRenorm.vmerge(s1,s1)
    s1 = CornerSpaceRenorm.hmerge(s1,s1)

    s2 = System(L,H,destroy(lb),J)
    s2 = CornerSpaceRenorm.vmerge(s2,s2)
    s2 = CornerSpaceRenorm.hmerge(s2,s2)

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
end
