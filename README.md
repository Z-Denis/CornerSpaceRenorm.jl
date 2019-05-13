# CornerSpaceRenorm.jl

[![Build Status](https://travis-ci.com/Z-Denis/CornerSpaceRenorm.jl.svg?token=XuYcpCDomapYmd2vHj9y&branch=master)](https://travis-ci.com/Z-Denis/CornerSpaceRenorm.jl)
[![codecov](https://codecov.io/gh/Z-Denis/CornerSpaceRenorm.jl/branch/master/graph/badge.svg?token=EwifsJO3ew)](https://codecov.io/gh/Z-Denis/CornerSpaceRenorm.jl)

In-development Julia package for performing corner space renormalization<sup id="a1">[1](#f1)</sup>.

## Principle

The corner space renormalization (CSR) procedure <sup id="a1">[1](#f1)</sup> consists in reducing the number of states required to represent a quantum mixed steady-state of a large lattice system by using previous knowledge about the steady-state of its subspaces. First, the steady-state <img src="https://latex.codecogs.com/gif.latex?\hat{\rho}" title="\hat{\rho}" /> of a reasonably sized lattice is computed, then two copies A and B of this lattice are merged and the CSR procedure is used to restrict the dimensionality of the doubled lattice Hilbert space <img src="https://latex.codecogs.com/gif.latex?\mathcal{H}_{A{\oplus}B}" title="\mathcal{H}_{A{\oplus}B}" /> from <img src="https://latex.codecogs.com/gif.latex?\lvert\mathcal{H}_{A}\rvert\times\lvert\mathcal{H}_{B}\rvert" title="\lvert\mathcal{H}_{A}\rvert\times\lvert\mathcal{H}_{B}\rvert" /> to some arbitrary corner dimension M, yielding the corner space <img src="https://latex.codecogs.com/gif.latex?\mathcal{H}_{A{\oplus}B}^{(M)}" title="\mathcal{H}_{A{\oplus}B}^{(M)}" />. Then the lattice Liouvillian is reconstructed on this compressed space. This process is then iterated by doubling the lattice over and over up to the desired size while monitoring the convergence of the compression method on some observables of interest.

### The merging procedure

From the steady-states <img src="https://latex.codecogs.com/gif.latex?\hat{\rho}_S&space;=&space;{\textstyle\sum_k}&space;p^{(S)}_k&space;\lvert\phi_k^{S}\rangle\langle\phi_k^{S}\rvert,\;&space;S\in\lbrace&space;A,B\rbrace" title="\hat{\rho}_S = {\textstyle\sum_k} p^{(S)}_k \lvert\phi_k^{S}\rangle\langle\phi_k^{S}\rvert,\; S\in\lbrace A,B\rbrace" /> of two compatible lattices, the Hilbert space of the doubled lattice <img src="https://latex.codecogs.com/gif.latex?\mathcal{H}_{A{\oplus}B}" title="\mathcal{H}_{A{\oplus}B}" /> can be given a basis <img src="https://latex.codecogs.com/gif.latex?\big\lbrace\lvert\phi_{r_k}^{A}\rangle\otimes\lvert\phi_{r_k^\prime}^{B}\rangle:&space;\;p_{k}^{(A)}p_{k}^{(B)}&space;\geq&space;p_{k&plus;1}^{(A)}p_{k&plus;1}^{(B)}&space;\big\rbrace_k" title="\big\lbrace\lvert\phi_{r_k}^{A}\rangle\otimes\lvert\phi_{r_k^\prime}^{B}\rangle: \;p_{k}^{(A)}p_{k}^{(B)} \geq p_{k+1}^{(A)}p_{k+1}^{(B)} \big\rbrace_k" /> of product states with decreasing joint probabilities. The corner space is then built as <img src="https://latex.codecogs.com/gif.latex?\mathcal{H}_{A{\oplus}B}^{(M)}&space;=&space;\mathrm{Span}\big(\big\lbrace\lvert\phi_{r_k}^{A}\rangle\otimes\lvert\phi_{r_k^\prime}^{B}\rangle&space;\big\rbrace_{k=1}^M\big)" title="\mathcal{H}_{A{\oplus}B}^{(M)} = \mathrm{Span}\big(\big\lbrace\lvert\phi_{r_k}^{A}\rangle\otimes\lvert\phi_{r_k^\prime}^{B}\rangle \big\rbrace_{k=1}^M\big)" /> and all system operators are mapped to this new space by means of a map

<img src="https://latex.codecogs.com/gif.latex?\mathcal{P}^{(M)}:&space;B_1(\mathcal{H}_{A})&space;\oplus&space;B_1(\mathcal{H}_{A})&space;\rightarrow&space;B_1\big(\mathcal{H}_{A{\oplus}B}^{(M)}\big),&space;\;\lbrace&space;\hat{O}_A,&space;\hat{O}_B\rbrace&space;\mapsto&space;\mathcal{P}^{(M)}(\hat{O}_A,&space;\hat{O}_B)" title="\mathcal{P}^{(M)}: B_1(\mathcal{H}_{A}) \oplus B_1(\mathcal{H}_{A}) \rightarrow B_1\big(\mathcal{H}_{A{\oplus}B}^{(M)}\big), \;\lbrace \hat{O}_A, \hat{O}_B\rbrace \mapsto \mathcal{P}^{(M)}(\hat{O}_A, \hat{O}_B)" />

defined as

<img src="https://latex.codecogs.com/gif.latex?\mathcal{P}^{(M)}(\hat{O}_A,&space;\hat{O}_B)&space;=&space;{\textstyle\sum_k}\langle\phi_{r_k}^{A}\rvert\hat{O}_A\lvert\phi_{r_{k^\prime}}^{A}\rangle\langle\phi_{r_k^\prime}^{B}\rvert\hat{O}_B\lvert\phi_{r_{k^\prime}^\prime}^{B}\rangle&space;\lvert\phi_{r_k}^{A},&space;\phi_{r_k^\prime}^{B}\rangle\langle\phi_{r_{k^\prime}}^{A},&space;\phi_{r_{k^\prime}^\prime}^{B}\rvert" title="\mathcal{P}^{(M)}(\hat{O}_A, \hat{O}_B) = {\textstyle\sum_k}\langle\phi_{r_k}^{A}\rvert\hat{O}_A\lvert\phi_{r_{k^\prime}}^{A}\rangle\langle\phi_{r_k^\prime}^{B}\rvert\hat{O}_B\lvert\phi_{r_{k^\prime}^\prime}^{B}\rangle \lvert\phi_{r_k}^{A}, \phi_{r_k^\prime}^{B}\rangle\langle\phi_{r_{k^\prime}}^{A}, \phi_{r_{k^\prime}^\prime}^{B}\rvert" />

For instance, the resulting Hamiltonian takes the form <img src="https://latex.codecogs.com/gif.latex?\hat{H}_{AB}^{(M)}&space;=&space;\mathcal{P}^{(M)}(\hat{H}_A,&space;\hat{1}_B)&space;&plus;&space;\mathcal{P}^{(M)}(\hat{1}_A,&space;\hat{H}_B)&space;&plus;&space;{\textstyle&space;\sum_{\langle&space;i\in&space;A;&space;j\in&space;B\rangle}}\mathcal{P}^{(M)}(\hat{O}_i^\dagger,&space;\hat{O}_j)&space;&plus;&space;\mathrm{H.c}" title="\hat{H}_{AB}^{(M)} = \mathcal{P}^{(M)}(\hat{H}_A, \hat{1}_B) + \mathcal{P}^{(M)}(\hat{1}_A, \hat{H}_B) + {\textstyle \sum_{\langle i\in A; j\in B\rangle}}\mathcal{P}^{(M)}(\hat{O}_i^\dagger, \hat{O}_j) + \mathrm{H.c}" />,

where <img src="https://latex.codecogs.com/gif.latex?\lbrace\hat{O}_j\rbrace_j" title="\lbrace\hat{O}_j\rbrace_j" /> are some coupling operators at the lattices A and B.

This package provides functions for applying this merging procedure in a multithreaded fashion for arbitrary N-dimensional lattices.

## Installation
From the Pkg REPL (prefix `]`), type:
```julia-repl
(v1.1) pkg> add https://github.com/Z-Denis/CornerSpaceRenorm.jl
```

## Lattices

####  SquareLattice (in deprecation process)

```julia
L = SquareLattice(nx,ny)
```
Construct an  `nx` times `ny` square lattice with open boundary conditions.
```julia
union(L1, L2)
L1 ∪ L2
```
Merge two square lattices along some compatible dimension.
```julia
vunion(L1, L2)
hunion(L1, L2)
```
Merge two square lattices either vertically or horizontally.

####  NdLattice

```julia
L = NdLattice((n1, n2, ... , nN); periodic=true/false)
```
Construct a ![alt text](https://latex.codecogs.com/gif.latex?\mathbb{Z}^N) lattice of size `n1 x  n2 x ... x nN` with periodic/open boundary conditions.
```julia
union(L1, L2, d)
```
Merge two N-dimensional lattices along the dimension `d`.

## Hamiltonian and dissipators

Hamiltonian and dissipators can be built from a lattice and some input operators.
```julia
hamiltonian(L, lH, t, lHt)
```
Generate a Hamiltonian for a lattice from some local Hamiltonian `lH` and a local hopping operator `lHt` with associated tunnelling rate `t` as <img src="https://latex.codecogs.com/gif.latex?{\textstyle&space;\sum_i}\hat{\texttt{lH}}_i&space;&plus;&space;{\textstyle&space;\sum_{\langle&space;i;j\rangle}}&space;(t\times\hat{\texttt{lHt}}_i\otimes\hat{\texttt{lHt}}_j^\dagger&space;&plus;&space;\mathrm{H.c.})" title="{\textstyle \sum_i}\hat{\texttt{lH}}_i + {\textstyle \sum_{\langle i;j\rangle}} (t\times\hat{\texttt{lHt}}_i\otimes\hat{\texttt{lHt}}_j^\dagger + \mathrm{H.c.})" />.
```julia
hamiltonian(L, lH, (t1,t2, ... ,tN), lHt)
```
Generate a Hamiltonian for a `NdLattice` in the same way but with anisotropic tunnelling rates.
```julia
dissipators(L, J)
```
Generate a vector of jump operators for the lattice from a vector of local jump operators.

## Systems

Corner merging and evolving functions take systems as arguments. Systems contain a lattice, a Hamiltonian, a set of jump operators, some other operators required during merging operations and a set of observables to be transformed into the corner subspace along with the Liouvillian.

####  SquareSystem (in deprecation process)
```julia
SquareSystem(L, H, (t, lHt), J)
SquareSystem(L, H, (t, lHt), J, obs)
```
Generate a system from a `SquareLattice`, a Hamiltonian, a tuple containing a tunnelling rate `t` and a local hopping operator `lHt` and a vector `J` of local jump operators. If passed, `obs` contains the observables to be transformed into the corner subspace along with the Liouvillian. Either a vector of `Dict` mapping some operator names (`String`) to some global operators or a vector of `Dict` mapping some operator names (`String`) to some local operators can be passed as argument.

####  NdSystem
```julia
NdSystem(L, H, t, lHt, J)
NdSystem(L, H, (t1,t2, ... ,tN), lHt, J)
NdSystem(L, H, t, J, lHt, obs)
NdSystem(L, H, (t1,t2, ... ,tN), J, lHt, obs)
```
Generate a system from a `NdLattice`, a Hamiltonian, either a single tunnelling rate `t` (isotropic system) or a tuple with as many rates as lattice dimensions (anisotropic system), and a local hopping operator `lHt` and a vector `J` of local jump operators. If passed, `obs` contains the observables to be transformed into the corner subspace along with the Liouvillian. Either a vector of `Dict` mapping some operator names (`String`) to some global operators or a vector of `Dict` mapping some operator names (`String`) to some local operators can be passed as argument.

## More generic Lattices

More generic lattices are accessible via the type `NdLattice{N}`. Its main constructor
interface `NdLattice(shape; periodic=true/false)` builds by default an N-dimensional cubic lattice, however it is possible to define more generic unit cells by directly
using the default constructor `NdLattice{N}(g, shape, Vint, Vext, pbc)`, where `N` is the number of dimensions, `g` is some arbitrary `LightGraphs.SimpleGraph` (unit cell), `shape` is the shape (`Tuple`) of the cubic lattice in which one wishes to embed the arbitrary unit cell,`Vint` is a tuple which i-th entry contains a vector of input vertices along the i-th dimension, `Vext` is a tuple which i-th entry contains a vector of output vertices along the i-th dimension and `pbc` is to be set to either `true` or `false` for either periodic or open boundary conditions. During the merging of two such lattices `L1` and `L2`, vertices belonging to `L1.Vext` are connected to `L2.Vint` (and reciprocally for `pbc=true`). The `shape` argument is used when `pbc=true`, for instance a one-dimensional periodic lattice `L=NdLattice((2,); periodic=true)` of size `2` (`L.shape[1] < 3`) has no supplementary edges with respect to its open boundary conditions counterpart (`NdLattice((2,); periodic=false)`) but once doubled `L2=union(L1,L1,1)` (`L2.shape[1] > 3`) a supplementary edge is added to complete the loop. This behaviour is taken care of similarly for generic lattices provided a proper shape was chosen at construction.

For instance, a Lieb lattice can be built as a two-dimensional `NdLattice` starting from a 3-vertex path graph and specifying custom input and output vertices for each merging direction:
```julia
using LightGraphs

g = PathGraph(3)
L = NdLattice{2}(g, (2,2), ([2],[2]), ([3],[1]), true)
```
Here `g` is a L shaped unit cell with vertices `1` (North), `2` (center) and `3`
(East). `(2, 2)` is a `2-Tuple` indicating that we embed this unit into a `2x2`
unit cell. `([2],[2])` (resp. `([3],[1])`) declares the input (resp. output)
vertices for the first and second dimension, i.e. if we append two such lattices
`L1` and `L2` along the dimension `1`, the merging would be operated by connecting the output vertex `3` of `L1` to the input vertex `2` of `L2`. In this way, if we would like to start the renormalization procedure from, e.g., a `2x2` Lieb lattice (12 vertices), we can just double `L` once along both directions:
```julia
L = union(L,L,1)
L = union(L,L,2)
```
We can then build a `ZnSystem` from it by simply using `NdSystem(L, H, t, lHt, J)` or any other `NdSystem` constructor.

## Corner operations

#### On SquareSystem instances (in deprecation process)
```julia
s1Us2 = merge(s1, s2)
```
Compression-free merging of two square systems along any compatible direction.

```julia
s1Us2 = vmerge(s1, s2)
```
Compression-free vertical merging of two square systems.

```julia
s1Us2 = hmerge(s1, s2)
```
Compression-free horizontal merging of two square systems.

```julia
s1Us2 = merge(s1, s2, ρ1, ρ2, M)
```
Projection into the corner space spanned by the `M` most probable product states of the merging of two square systems along any compatible direction.

```julia
s1Us2 = vmerge(s1, s2, ρ1, ρ2, M)
```
Projection into the corner space spanned by the `M` most probable product states of the merging of two square systems along the vertical direction.

```julia
s1Us2 = hmerge(s1, s2, ρ1, ρ2, M)
```
Projection into the corner space spanned by the `M` most probable product states of the merging of two square systems along the horizontal direction.

#### On NdSystem instances

```julia
s1Us2 = merge(s1, s2, d, ρ1, ρ2, M)
```
Projection into the corner space spanned by the `M` most probable product states of the merging of two N-dimensional systems along the direction `d`.

## Steady state evaluation

```julia
ρ = steadystate_bicg(H, J, l; log=false, kwargs...)
ρ, log = steadystate_bicg(H, J, l; log=true, kwargs...)
ρ = steadystate_bicg(s, l; log=false, kwargs...)
ρ, log = steadystate_bicg(s, l; log=true, kwargs...)
```
Evaluates the steady state by solving iteratively the linear system <img src="https://latex.codecogs.com/gif.latex?\mathcal{L}\hat{\rho}&space;=&space;0" title="\mathcal{L}\hat{\rho} = 0" /> via the stabilized biconjugate gradient method with `l` `GMRES` steps. The first line of the Liouvillian is overwritten to enforce a non-trivial trace one solution, this approximation yields an error of the order of the inverse of the square of the size of the Hilbert space.

```julia
ρ = steadystate_bicg_LtL(H, J, l; log=false, kwargs...)
ρ, log = steadystate_bicg_LtL(H, J, l; log=true, kwargs...)
ρ = steadystate_bicg_LtL(s, l; log=false, kwargs...)
ρ, log = steadystate_bicg_LtL(s, l; log=true, kwargs...)
```
Evaluates the steady state by solving iteratively the linear system <img src="https://latex.codecogs.com/gif.latex?\langle\mathcal{L},\mathcal{L}\hat{\rho}\rangle&space;&plus;&space;\langle\mathrm{Tr},\hat{\rho}\rangle&space;\mathrm{Tr}&space;=&space;\mathrm{Tr}" title="\langle\mathcal{L},\mathcal{L}\hat{\rho}\rangle + \langle\mathrm{Tr},\hat{\rho}\rangle \mathrm{Tr} = \mathrm{Tr}" /> via the stabilized biconjugate gradient method with `l` `GMRES` steps. No approximation of the Liovillian is made in order to enfore the trace one but in practice convergence is slower and poorer.

## Plotting

Methods for `Lattice` and `AbstractSystem` are added to function `GraphPlot.gplot` of the package [GraphPlot.jl](https://github.com/JuliaGraphs/GraphPlot.jl).
```julia
GraphPlot.gplot(L; kwargs...)
GraphPlot.gplot(L, locs_x_in, locs_y_in; kwargs...)
GraphPlot.gplot(s; kwargs...)
GraphPlot.gplot(s, locs_x_in, locs_y_in; kwargs...)
plot_system(s; kwargs...)
plot_system(s, locs_x_in, locs_y_in; kwargs...)
```

## Examples

#### 16-site dissipative quantum transverse Ising model

```julia
using CornerSpaceRenorm, QuantumOptics

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

# For simplicity we target an easy (low entropy) regime of parameters
g = 1.          # Local energies
gamma = 1. # Dissipation rate
V = 2.          # tunneling rate

# 1D lattice with periodic boundary conditions. Dimension is guessed from the
# length of the shape tuple.
L = NdLattice((8,); periodic=true)
H = hamiltonian(L, g/2 * sx, V/4, sz)
J = dissipators(L, [sqrt(2gamma) * sm])

# Generate a system from a lattice, a Hamiltonian,
# a local tunnelling operator and jump operators
s = NdSystem(L, H, V/4., sz, J, lobs)
# Compute the steady state (by brute-force integration)
ρ = steadystate.master(s)[2][end]
# Merge two systems into some corner subspace spanned by 500 kets
s2 = merge(s,s,1,ρ,ρ,500) # 16 x 1 lattice
# Find the steady-state density matrix by minimizing Lρ via the BiCGStab(l=4)
# iterative method.
ρ2 = steadystate_bicg(s2, 4; tol=1e-4, verbose=true)
hermitianize!(ρ2)

# Compute some observables averaged over the entire lattice
σx_mean = real(sum([expect(s2.observables[i]["sigmax"],ρ2) for i in 1:16]) / 16)
σy_mean = real(sum([expect(s2.observables[i]["sigmay"],ρ2) for i in 1:16]) / 16)
σz_mean = real(sum([expect(s2.observables[i]["sigmaz"],ρ2) for i in 1:16]) / 16)

# Print some results
println("Purity = ",real(tr(ρ2*ρ2)),"\n#states = ", real(exp(entropy_vn(ρ2))))
println("<sx> = $σx_mean\n<sy> = $σy_mean\n<sz> = $σz_mean")
```
```julia-repl
julia> Purity = 0.9741461232157818
julia> #states = 1.1184964234832833
julia> <sx> = 0.23599186171061198
julia> <sy> = 0.062087577378825566
julia> <sz> = -0.966768810606053
```
Increasing the corner dimension from `500` to `700` yields
```julia-repl
julia> Purity = Purity = 0.9741118335231511
julia> #states = 1.1185224207488438
julia> <sx> = 0.2359865263666315
julia> <sy> = 0.06208640000017497
julia> <sz> = -0.9667527388195073
```
The corner has reached convergence up to the chosen tolerance.

#### 8x8 hardcore-boson driven-dissipative lattice

```julia
using CornerSpaceRenorm, QuantumOptics

# Local basis
b = FockBasis(1)

# Build some operators in the local basis
N = number(b)
a = destroy(b)

# Build a set of local observables in the local basis
lobs = Dict("N"=>N)

κ = 1     # Dissipation rate
F = 1e-1κ # Driving strength
t = 1.    # tunneling rate
Δ = -1t   # Detuning

# 2x4 Bose-Hubbard lattice with periodic boundary conditions.
# Dimension is guessed from the length of the shape tuple.
L = NdLattice((2,4); periodic=true)
# Local Hamiltonian
lH = -Δ*N + F*(a + dagger(a))
H = hamiltonian(L, lH, -t/2., a)
J = dissipators(L, [sqrt(κ) * a])

# Generate a system from a lattice, a Hamiltonian,
# a local tunnelling operator and jump operators
s = NdSystem(L, H, -t/2., a, J, lobs)
# Compute the steady state (by brute-force integration)
ρ = steadystate.master(s)[2][end]
# Merge two systems into some corner subspace spanned by 500 kets
s2 = merge(s,s,1,ρ,ρ,500) # 4x4-site lattice
# Find the steady-state density matrix by minimizing Lρ via the BiCGStab(l=4)
# iterative method.
ρ2 = steadystate_bicg(s2, 6; tol=1e-4, verbose=true)
hermitianize!(ρ2)
# Repeat the process
s3 = merge(s2,s2,1,ρ2,ρ2,500) # 8x4-site lattice
ρ3 = steadystate_bicg(s3, 6; tol=1e-4, verbose=true)
hermitianize!(ρ3)
# etc.
s4 = merge(s3,s3,2,ρ3,ρ3,700) # 8x8-site lattice
ρ4 = steadystate_bicg(s4, 6; tol=1e-4, verbose=true)
hermitianize!(ρ4)

mean_pop(s,ρ) = real(sum([expect(s.observables[i]["N"],ρ) for i in 1:nv(s.lattice)]) / nv(s.lattice))

# Print some results
println("\nIteration 1")
println("Purity = ",real(tr(ρ*ρ)),"\n#states = ", real(exp(entropy_vn(ρ))))
println("<N> = ",mean_pop(s,ρ))

println("\nIteration 2")
println("Purity = ",real(tr(ρ2*ρ2)),"\n#states = ", real(exp(entropy_vn(ρ2))))
println("<N> = ",mean_pop(s2,ρ2))

println("\nIteration 3")
println("Purity = ",real(tr(ρ3*ρ3)),"\n#states = ", real(exp(entropy_vn(ρ3))))
println("<N> = ",mean_pop(s3,ρ3))

println("\nIteration 4")
println("Purity = ",real(tr(ρ4*ρ4)),"\n#states = ", real(exp(entropy_vn(ρ4))))
println("<N> = ",mean_pop(s4,ρ4))
```
```julia-repl
julia>
Iteration 1
Purity = 0.9856262773126672
#states = 1.0580681393505547
<N> = 0.02135711148025299

Iteration 2
Purity = 0.9380319564966497
#states = 1.2614296231492617
<N> = 0.011410774377058204

Iteration 3
Purity = 0.619858327346652
#states = 3.8137560169203337
<N> = 0.022882387269542788

Iteration 4
Purity = 0.03660696744651538
#states = 130.393737931304
<N> = 0.07190099978628699
```

## References

<b id="f1">[1]</b> S. Finazzi, A. Le Boité, F. Storme, A. Baksic, and C. Ciuti, Corner-Space Renormalization Method for Driven-Dissipative Two-Dimensional Correlated Systems, [Phys. Rev. Lett. 115, 080604 (2015)](https://doi.org/10.1103/PhysRevLett.115.080604) [↩](#a1)
