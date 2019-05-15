using CornerSpaceRenorm, QuantumOptics
using LightGraphs, GraphPlot
using LinearAlgebra, IterativeSolvers, Random
using Test, InteractiveUtils

@testset "CornerSpaceRenorm" begin
    @testset "Operators" begin
        lb = FockBasis(1)
        lH = randoperator(lb)
        tH = randoperator(lb)

        for pbc in [true, false]
            L = NdLattice((3,2,1);periodic=false)
            H1 = hamiltonian(L,lH,(1.,1.,1.),tH)
            H2 = hamiltonian(L,lH,(1.,tH))
            H3 = hamiltonian(L,lH,1.,tH)
            @test maximum(abs2.(H1.data .- H2.data)) ≈ 0. atol=eps(Float64)
            @test maximum(abs2.(H1.data .- H3.data)) ≈ 0. atol=eps(Float64)
        end
    end;

    @testset "Types" begin
        @testset "NdLattice" begin
            L = NdLattice((2,2))

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

            function compare_lattices(L1,L2)
                @test typeof(L1) == typeof(L2)
                @test all([getfield(L1,f) == getfield(L2,f) for f in fieldnames(typeof(L1))])
            end

            # Test lattice merging
            L1 = NdLattice((2,3))
            L2 = NdLattice((2,1))
            L3_12 = union(L1,L2,2)
            L3_21 = union(L2,L1,2)
            compare_lattices(L3_12,L3_21)
            @test_throws DimensionMismatch union(L1,L2,1)

            L1 = NdLattice((3,2))
            L2 = NdLattice((1,2))
            L3_12 = union(L1,L2,1)
            L3_21 = union(L2,L1,1)
            @test typeof(L3_12) == typeof(L3_21)
            # TO DO: test boundaries
            @test LightGraphs.Experimental.could_have_isomorph(L3_12.L, L3_21.L)
            @test L3_12.shape[1] == L1.shape[1] + L2.shape[1]
            compare_lattices(L3_12, union(L1,L2,1))
            @test_throws DimensionMismatch union(L1,L2,2)

            Lpbc = NdLattice((2,2);periodic=true)
            Lobc = NdLattice((2,2);periodic=false)
            @test_throws AssertionError union(Lpbc,Lobc,1)
            @test_throws AssertionError union(Lpbc,Lobc,2)

            L = NdLattice((2,);periodic=true)
            LuL = union(L,L,1)
            L2 = NdLattice((4,);periodic=true)
            compare_lattices(LuL,L2)

            L = NdLattice((2,2);periodic=true)
            LuL = union(L,L,1)
            L2 = NdLattice((4,2);periodic=true)
            @test typeof(LuL) == typeof(L2)
            @test LightGraphs.Experimental.could_have_isomorph(LuL.L, L2.L)

            L = NdLattice((2,2);periodic=true)
            LuL = union(L,L,2)
            L2 = NdLattice((2,4);periodic=true)
            @test typeof(LuL) == typeof(L2)
            compare_lattices(LuL,L2)

            L = NdLattice((1,);periodic=true)
            for i in 1:3
                L = union(L,L,1)
            end
            compare_lattices(L, NdLattice((8,); periodic=true))

            # Dummy tests on plotting
            x = gplot(L)
            y = gplot(L.L)
            @test typeof(x) == typeof(y)
            x = gplot(L, Float64.(collect(1:nv(L))), Float64.(collect(1:nv(L))))
            y = gplot(L.L, Float64.(collect(1:nv(L))), Float64.(collect(1:nv(L))))
            @test all([typeof(getfield(x,f)) == typeof(getfield(y,f)) for f in fieldnames(typeof(x))])
        end;

        @testset "NdSystem" begin
            # Local basis
            b_spin = SpinBasis(1//2)

            # Build some operators in the local basis
            sm = sigmam(b_spin)
            sp = sigmap(b_spin)
            sz = sigmaz(b_spin)
            sx = sigmax(b_spin)
            sy = sigmay(b_spin)

            # Build a set of local observables in the local basis
            lobs = Dict("sigmax"=>sx,"sigmay"=>sy,"sigmaz"=>sz)

            g = 1.
            gamma = 1.
            V = 2/2.

            # Test assertions in constructors
            L = NdLattice((4,))
            H = hamiltonian(L, g/2 * sx, V/4, sz)
            J = dissipators(L, [sqrt(2gamma) * sm])
            # Incompatible lattice
            L_err = NdLattice((5,))
            @test_throws AssertionError NdSystem(L_err, H, (V/4,), sz, J, lobs)
            # Number of tunnelling rates incompatible with number of dimensions
            @test_throws AssertionError NdSystem(L, H, (V/4,V/4), sz, J, lobs)
            # Incompatible operators
            @test_throws MethodError NdSystem(L, H, (V/4,), destroy(FockBasis(1)), J, lobs)

            # Test constructor with missing observables
            s = NdSystem(L, H, (V/4,), sz, J)
            @test eltype(s.observables) == Dict{String,typeof(H)}

            # Broken tests on displaying lattices
            @test_skip gplot(L) == gplot(L.L)
            @test_skip plot_system(s) == gplot(s.lattice.L)

            # Broken: cannot construct a system with 1 site.

            function compare_systems(s1::AbstractSystem,s2::AbstractSystem)
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
            # Simpler constructor for homogeneous coupling
            L = NdLattice((2,1,1,1))
            H = hamiltonian(L, g/2 * sx, V/4, sz)
            J = dissipators(L, [sqrt(2gamma) * sm])
            s1 = NdSystem(L, H, V/4, sz, J, lobs)
            s2 = NdSystem(L, H, Tuple([V/4 for i in 1:4]), sz, J, lobs)
            compare_systems(s1,s2)
            s1 = NdSystem(L, H, V/4, sz, J)
            s2 = NdSystem(L, H, Tuple([V/4 for i in 1:4]), sz, J)
            compare_systems(s1,s2)

            # Test construction and merging
            # 1D
            for pbc in [true, false]
                L = NdLattice((4,); periodic=pbc)
                H = hamiltonian(L, g/2 * sx, V/4, sz)
                J = dissipators(L, [sqrt(2gamma) * sm])
                s = NdSystem(L, H, (V/4,), sz, J, lobs)
                ρ = steadystate.master(s)[2][end]
                sUs = merge(s,s,1,ρ,ρ,256)

                L2 = NdLattice((8,); periodic=pbc)
                H2 = hamiltonian(L2, g/2 * sx, V/4, sz)
                J2 = dissipators(L2, [sqrt(2gamma) * sm])
                s2 = NdSystem(L2, H2, (V/4,), sz, J2, lobs)

                @test maximum(abs2.(eigvals(sUs.H.data) .- eigvals(Matrix(s2.H.data)))) ≈ 0. atol=eps(Float64)
                @test maximum([maximum(abs2.(eigvals(sUs.J[i].data'sUs.J[i].data) .- eigvals(Matrix(s2.J[i].data's2.J[i].data)))) for i in 1:length(sUs.J)]) ≈ 0. atol=eps(Float64)

                d = 1
                @test maximum([maximum(abs2.(eigvals(sUs.Htint[d][i].data) .- eigvals(Matrix(s2.Htint[d][i].data)))) for i in 1:length(sUs.Htint[d])]) ≈ 0. atol=eps(Float64)
                @test maximum([maximum(abs2.(eigvals(sUs.Htext[d][i].data) .- eigvals(Matrix(s2.Htext[d][i].data)))) for i in 1:length(sUs.Htext[d])]) ≈ 0. atol=eps(Float64)
            end

            # 2D
            for pbc in [true, false]
                L = NdLattice((2,2); periodic=pbc)
                H = hamiltonian(L, g/2 * sx, V/4, sz)
                J = dissipators(L, [sqrt(2gamma) * sm])
                s = NdSystem(L, H, (V/4, V/4), sz, J, lobs)
                ρ = steadystate.master(s)[2][end]
                sUs = merge(s,s,2,ρ,ρ,256)

                L2 = NdLattice((2,4); periodic=pbc)
                H2 = hamiltonian(L2, g/2 * sx, V/4, sz)
                J2 = dissipators(L2, [sqrt(2gamma) * sm])
                s2 = NdSystem(L2, H2, (V/4, V/4), sz, J2, lobs)

                @test maximum(abs2.(eigvals(sUs.H.data) .- eigvals(Matrix(s2.H.data)))) ≈ 0. atol=eps(Float64)
                @test maximum([maximum(abs2.(eigvals(sUs.J[i].data'sUs.J[i].data) .- eigvals(Matrix(s2.J[i].data's2.J[i].data)))) for i in 1:length(sUs.J)]) ≈ 0. atol=eps(Float64)

                for d in 1:2
                    @test maximum([maximum(abs2.(eigvals(sUs.Htint[d][i].data) .- eigvals(Matrix(s2.Htint[d][i].data)))) for i in 1:length(sUs.Htint[d])]) ≈ 0. atol=eps(Float64)
                    @test maximum([maximum(abs2.(eigvals(sUs.Htext[d][i].data) .- eigvals(Matrix(s2.Htext[d][i].data)))) for i in 1:length(sUs.Htext[d])]) ≈ 0. atol=eps(Float64)
                end
            end
            # TO DO: Test construction and merging of generic lattices

            # Dummy tests on plotting
            L = NdLattice((4,))
            s = NdSystem(L, H, (V/4,), sz, J)
            x = gplot(s)
            y = gplot(L.L)
            @test typeof(x) == typeof(y)
            x = plot_system(s)
            @test typeof(x) == typeof(y)
            x = gplot(s, Float64.(collect(1:nv(L))), Float64.(collect(1:nv(L))))
            y = gplot(L.L, Float64.(collect(1:nv(L))), Float64.(collect(1:nv(L))))
            @test all([typeof(getfield(x,f)) == typeof(getfield(y,f)) for f in fieldnames(typeof(x))])
            x = plot_system(s, Float64.(collect(1:nv(L))), Float64.(collect(1:nv(L))))
            @test all([typeof(getfield(x,f)) == typeof(getfield(y,f)) for f in fieldnames(typeof(x))])
        end;

        @testset "CornerBasis" begin
            # Test constructors
            b1 = SpinBasis(1//2)
            b2 = CompositeBasis([b1 for i in 1:10])
            bC = CornerBasis(b1, b2, 10)
            @test bC == CornerBasis([10],10,[[2], fill(2,10)])
            bC2 = CornerBasis(bC, bC, 10)
            @test bC2 == CornerBasis([10],10,vcat(bC.local_shapes, bC.local_shapes))

            # Test Base.== overload
            bC = CornerBasis(10)
            @test bC == deepcopy(bC)
            @test bC !== SpinBasis(1//2)
            @test bC !== FockBasis(1)
            @test bC !== GenericBasis(1)
            bC2 = CornerBasis([11],bC.M,bC.local_shapes)
            @test bC !== bC2
            bC2 = CornerBasis(bC.shape,bC.M+1,bC.local_shapes)
            @test bC !== bC2
            bC2 = CornerBasis(bC.shape,bC.M,vcat(bC.local_shapes, [[1]]))
            @test bC !== bC2
        end;
    end;

    @testset "Corner algorithm" begin
        @testset "Core" begin
            # Test orthonormalisation
            N = 10
            U = diagm([i=>ones(ComplexF64,N-i) for i in 0:N]...)
            @test maximum(abs2.(CornerSpaceRenorm.cgs(U) .- Diagonal(ones(ComplexF64,N)))) ≈ 0. atol=eps(Float64)
            @test maximum(abs2.(CornerSpaceRenorm.mgs(U) .- Diagonal(ones(ComplexF64,N)))) ≈ 0. atol=eps(Float64)
            V = copy(U); CornerSpaceRenorm.cgs!(V)
            @test maximum(abs2.(V .- Diagonal(ones(ComplexF64,N)))) ≈ 0. atol=eps(Float64)
            V = copy(U); CornerSpaceRenorm.mgs!(V)
            @test maximum(abs2.(V .- Diagonal(ones(ComplexF64,N)))) ≈ 0. atol=eps(Float64)

            # Test maxk
            x = randperm(1000)
            idcs, y = CornerSpaceRenorm.maxk(x,100)
            # Test that indices match sorted values
            @test all([x[idcs[i]] == y[i] for i in 1:100])
            # Test that k sorted values are the largest
            @test all([val == y[i] for (i,val) in enumerate(1000:-1:901)])

            # Test max_prod_pairs
            x1 = 1:50; x2 = 1:50
            idcs, pairs = CornerSpaceRenorm.max_prod_pairs(x1,x2,10)
            @test all([x1[idx[1]] * x2[idx[2]] for idx in idcs] .== pairs)
            @test all(sort([x1[i] * x2[j] for i in 1:50, j in 1:50][:],rev=true)[1:10] .== pairs)

            # Test hermitianize
            b = GenericBasis(50)
            Op = randoperator(b)
            OpH = hermitianize(Op)
            @test ishermitian(OpH.data)
            @test maximum(abs2.(OpH.data .- (Op + dagger(Op)).data/2.)) ≈ 0. atol=eps(Float64)
            Op = randoperator(b);
            OpH = copy(Op);
            hermitianize!(OpH)
            @test ishermitian(OpH.data)
            @test maximum(abs2.(OpH.data .- (Op + dagger(Op)).data/2.)) ≈ 0. atol=eps(Float64)
        end;

        @testset "Merging" begin
            @testset "NdSystem" begin
                # Local basis
                b_spin = SpinBasis(1//2)

                # Build some operators in the local basis
                sm = sigmam(b_spin)
                sp = sigmap(b_spin)
                sz = sigmaz(b_spin)
                sx = sigmax(b_spin)
                sy = sigmay(b_spin)

                # Build a set of local observables in the local basis
                lobs = Dict("sigmax"=>sx,"sigmay"=>sy,"sigmaz"=>sz)

                g = 1.
                gamma = 1.
                V = 2/2.

                L1 = NdLattice((2,1))
                H1 = hamiltonian(L1, g/2 * sx, V/4, sz)
                J1 = dissipators(L1, [sqrt(2gamma) * sm])
                s1 = NdSystem(L1, H1, (V/4,V/4), sz, J1, lobs)
                ρ1 = steadystate.master(s1.H,s1.J)[2][end]
                L2 = NdLattice((1,2))
                H2 = hamiltonian(L2, g/2 * sx, V/4, sz)
                J2 = dissipators(L2, [sqrt(2gamma) * sm])
                s2 = NdSystem(L2, H2, (V/4,V/4), sz, J2, lobs)
                ρ2 = steadystate.master(s2.H,s2.J)[2][end]

                @test_throws DimensionMismatch merge(s1, s2, 1, ρ1, ρ2, 2)
                @test_throws DimensionMismatch merge(s1, s2, 2, ρ1, ρ2, 2)

                L1_pbc = NdLattice((2,1);periodic=true)
                H1_pbc = hamiltonian(L1_pbc, g/2 * sx, V/4, sz)
                J1_pbc = dissipators(L1_pbc, [sqrt(2gamma) * sm])
                s1_pbc = NdSystem(L1_pbc, H1_pbc, (V/4,V/4), sz, J1_pbc, lobs)
                ρ1_pbc = steadystate.master(s1_pbc.H,s1_pbc.J)[2][end]

                @test_throws AssertionError merge(s1,s1_pbc,1,ρ1,ρ1_pbc,2)

                # Full space corner vs exact system
                # 1D
                for pbc in [true, false]
                    L1 = NdLattice((2,))
                    H1 = hamiltonian(L1, g/2 * sx, V/4, sz)
                    J1 = dissipators(L1, [sqrt(2gamma) * sm])
                    s1 = NdSystem(L1, H1, (V/4,), sz, J1, lobs)
                    ρ1 = steadystate.master(s1.H,s1.J)[2][end]
                    s1 = merge(s1,s1,1,ρ1,ρ1,16)
                    ρ1 = steadystate.master(s1.H,s1.J;tol=1e-5)[2][end]

                    L2 = NdLattice((4,))
                    H2 = hamiltonian(L2, g/2 * sx, V/4, sz)
                    J2 = dissipators(L2, [sqrt(2gamma) * sm])
                    s2 = NdSystem(L2, H2, (V/4,), sz, J2, lobs)
                    ρ2 = steadystate.master(s2.H,s2.J;tol=1e-5)[2][end]

                    @test maximum(abs2.(eigvals(Matrix(s1.H.data)) .- eigvals(Matrix(s2.H.data)))) ≈ 0. atol=eps(Float64)
                    @test maximum([maximum(abs2.(eigvals(Matrix(s1.J[i].data's1.J[i].data)) .- eigvals(Matrix(s2.J[i].data's2.J[i].data)))) for i in 1:length(s1.J)]) ≈ 0. atol=eps(Float64)
                    @test maximum(abs2.(eigvals(ρ1.data) .- eigvals(ρ2.data))) ≈ 0. atol=1e-5
                end

                # 2D
                for pbc in [true, false]
                    L1 = NdLattice((2,1))
                    H1 = hamiltonian(L1, g/2 * sx, V/4, sz)
                    J1 = dissipators(L1, [sqrt(2gamma) * sm])
                    s1 = NdSystem(L1, H1, (V/4,V/4), sz, J1, lobs)
                    ρ1 = steadystate.master(s1.H,s1.J)[2][end]
                    s1 = merge(s1,s1,1,ρ1,ρ1,16)
                    ρ1 = steadystate.master(s1.H,s1.J;tol=1e-5)[2][end]

                    L2 = NdLattice((4,1))
                    H2 = hamiltonian(L2, g/2 * sx, V/4, sz)
                    J2 = dissipators(L2, [sqrt(2gamma) * sm])
                    s2 = NdSystem(L2, H2, (V/4,V/4), sz, J2, lobs)
                    ρ2 = steadystate.master(s2.H,s2.J;tol=1e-5)[2][end]

                    @test maximum(abs2.(eigvals(Matrix(s1.H.data)) .- eigvals(Matrix(s2.H.data)))) ≈ 0. atol=eps(Float64)
                    @test maximum([maximum(abs2.(eigvals(Matrix(s1.J[i].data's1.J[i].data)) .- eigvals(Matrix(s2.J[i].data's2.J[i].data)))) for i in 1:length(s1.J)]) ≈ 0. atol=eps(Float64)
                    @test maximum(abs2.(eigvals(ρ1.data) .- eigvals(ρ2.data))) ≈ 0. atol=1e-5
                end

                # Merging different lattices
                Ls = NdLattice((3,))
                Hs = hamiltonian(Ls, g/2 * sx, V/4, sz)
                Js = dissipators(Ls, [sqrt(2gamma) * sm])
                ss = NdSystem(Ls, Hs, V/4, sz, Js, lobs)
                ρs = steadystate.master(ss)[2][end]

                t = 1.
                Δ = 1.
                κ = 1.
                F = 5e-1 * κ

                lb = FockBasis(1)
                Lb = NdLattice((3,),periodic=false)
                Hb = hamiltonian(Lb, -Δ * number(lb) + F * (create(lb) + destroy(lb)), t/2., destroy(lb))
                Jb = dissipators(Lb,[sqrt(κ) * destroy(lb)])
                sb = NdSystem(Lb,Hb,(t/2.,),destroy(lb),Jb)
                ρb = steadystate.master(sb)[2][end]

                @test_logs (:info,"The two input systems do not have equal tunnelling rates along all directions.") s_sb = merge(ss,sb,1,ρs,ρb,10)
            end;
        end;
    end;

    @testset "steadystate" begin
        # Local basis
        b_spin = SpinBasis(1//2)

        # Build some operators in the local basis
        sm = sigmam(b_spin)
        sp = sigmap(b_spin)
        sz = sigmaz(b_spin)
        sx = sigmax(b_spin)
        sy = sigmay(b_spin)

        # Build a set of local observables in the local basis
        lobs = Dict("sigmax"=>sx,"sigmay"=>sy,"sigmaz"=>sz)

        g = 1.
        gamma = 1.
        V = 2/2.

        L1 = NdLattice((6,))
        H1 = hamiltonian(L1, g/2 * sx, V/4, sz)
        J1 = dissipators(L1, [sqrt(2gamma) * sm])
        s1 = NdSystem(L1, H1, (V/4,), sz, J1, lobs)
        ρ1 = steadystate.master(s1.H,s1.J; tol=1e-5)[2][end]
        ρ2 = steadystate_bicg(s1)
        @test fidelity(ρ1, ρ2) ≈ 1. atol=1e-5
        ρ3, log = steadystate_bicg(s1; log=true)
        @test typeof(log) == ConvergenceHistory{true,Nothing}
        ρ4 = CornerSpaceRenorm.steadystate_bicg_LtL(s1)
        @test fidelity(ρ1, ρ4) ≈ 1. atol=1e-5
        ρ, log = steadystate_bicg(s1; log=true)
        @test typeof(log) == ConvergenceHistory{true,Nothing}
        ρ, log = CornerSpaceRenorm.steadystate_bicg_LtL(s1; log=true)
        @test typeof(log) == ConvergenceHistory{true,Nothing}
    end;

end;
