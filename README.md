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

<img src="https://latex.codecogs.com/gif.latex?\mathcal{P}^{(M)}(\hat{O}_A,&space;\hat{O}_B)&space;=&space;{\textstyle\sum_{k,k^\prime=1}^M}\langle\phi_{r_k}^{A}\rvert\hat{O}_A\lvert\phi_{r_{k^\prime}}^{A}\rangle\langle\phi_{r_k^\prime}^{B}\rvert\hat{O}_B\lvert\phi_{r_{k^\prime}^\prime}^{B}\rangle&space;\lvert\phi_{r_k}^{A},&space;\phi_{r_k^\prime}^{B}\rangle\langle\phi_{r_{k^\prime}}^{A},&space;\phi_{r_{k^\prime}^\prime}^{B}\rvert" title="\mathcal{P}^{(M)}(\hat{O}_A, \hat{O}_B) = {\textstyle\sum_{k,k^\prime=1}^M}\langle\phi_{r_k}^{A}\rvert\hat{O}_A\lvert\phi_{r_{k^\prime}}^{A}\rangle\langle\phi_{r_k^\prime}^{B}\rvert\hat{O}_B\lvert\phi_{r_{k^\prime}^\prime}^{B}\rangle \lvert\phi_{r_k}^{A}, \phi_{r_k^\prime}^{B}\rangle\langle\phi_{r_{k^\prime}}^{A}, \phi_{r_{k^\prime}^\prime}^{B}\rvert" />

For instance, the resulting Hamiltonian takes the form <img src="https://latex.codecogs.com/gif.latex?\hat{H}_{AB}^{(M)}&space;=&space;\mathcal{P}^{(M)}(\hat{H}_A,&space;\hat{1}_B)&space;&plus;&space;\mathcal{P}^{(M)}(\hat{1}_A,&space;\hat{H}_B)&space;&plus;&space;{\textstyle&space;\sum_{\langle&space;i\in&space;A;&space;j\in&space;B\rangle}}\mathcal{P}^{(M)}(\hat{O}_i^\dagger,&space;\hat{O}_j)&space;&plus;&space;\mathrm{H.c}" title="\hat{H}_{AB}^{(M)} = \mathcal{P}^{(M)}(\hat{H}_A, \hat{1}_B) + \mathcal{P}^{(M)}(\hat{1}_A, \hat{H}_B) + {\textstyle \sum_{\langle i\in A; j\in B\rangle}}\mathcal{P}^{(M)}(\hat{O}_i^\dagger, \hat{O}_j) + \mathrm{H.c}" />,

where <img src="https://latex.codecogs.com/gif.latex?\lbrace\hat{O}_j\rbrace_j" title="\lbrace\hat{O}_j\rbrace_j" /> are some coupling operators at the lattices A and B.

This package provides functions for applying this merging procedure in a multithreaded fashion for arbitrary N-dimensional lattices.

## Installation
From the Pkg REPL (prefix `]`), type:
```julia-repl
(v1.1) pkg> add https://github.com/Z-Denis/CornerSpaceRenorm.jl
```

## Lattices

#### NdLattice

```julia
L = NdLattice((n1, n2, ... , nN); periodic=true/false)
```
Construct a ![alt text](https://latex.codecogs.com/gif.latex?\mathbb{Z}^N) lattice of size `n1 x  n2 x ... x nN` with periodic/open boundary conditions.
```julia
union(L1, L2, d)
```
Merge two N-dimensional lattices along the dimension `d`.
```julia
pbc_from_obc(L)
```
Get a periodic copy of `L`.
```julia
obc_from_pbc(L)
```
Get an open copy of `L`.

#### More generic lattices

The corner space renormalization method was successfully applied to non-cubic lattices such as Lieb lattices<sup id="a2">[2](#f2)</sup> . Even though the above-presented `NdLattice` constructor builds by default N-dimensional cubic lattices, it is possible to define generic unit cells by calling directly the default constructor `NdLattice{N}(g, shape, Vint, Vext, pbc)`, where `N` is the number of dimensions, `g` is some arbitrary `LightGraphs.SimpleGraph` (unit cell), `shape` is the shape (`Tuple`) of the cubic lattice in which one wishes to embed the arbitrary unit cell,`Vint` is a tuple which i-th entry contains a vector of input vertices along the i-th dimension, `Vext` is a tuple which i-th entry contains a vector of output vertices along the i-th dimension and `pbc` is to be set to either `true` or `false` for either periodic or open boundary conditions. During the merging of two such lattices `L1` and `L2`, vertices belonging to `L1.Vext` are connected to `L2.Vint` (and reciprocally for `pbc=true`). The `shape` argument is used when `pbc=true`, for instance a one-dimensional periodic lattice `L=NdLattice((2,); periodic=true)` of size `2` (`L.shape[1] < 3`) has no supplementary edges with respect to its open boundary conditions counterpart (`NdLattice((2,); periodic=false)`) but once doubled `L2=union(L1,L1,1)` (`L2.shape[1] > 3`) a supplementary edge is added to complete the loop. This behavior is taken care of similarly for generic lattices provided a proper shape was chosen at construction.

For instance, let us consider a Lieb lattice which unit cell is defined as follows:

<p align="center">
  <img src="https://github.com/Z-Denis/CornerSpaceRenorm.jl/blob/development/fig_Lieb.png" align="center" height="200" width="200" ></a>
</p>

such a lattice can be built as a two-dimensional `NdLattice` starting from a 3-vertex path graph and specifying custom input and output vertices for each merging direction:

```julia
using LightGraphs

D = 2                # #dimensions
a = 1; b = 2; c = 3; # vertex names
shape = (2, 2)       # shape of the unit cell
pbc = true           # Boundary conditions
Vin_1  = [b]         # Input vertices in direction d=1
Vout_1 = [c]         # Output vertices in direction d=1
Vin_2  = [b]         # Input vertices in direction d=2
Vout_2 = [a]         # Output vertices in direction d=2
unit_cell = PathGraph(3)

L = NdLattice{D}(unit_cell, shape, (Vin_1, Vin_2), (Vout_1, Vout_2), pbc)
```

In this example, one declares the input (resp. output) vertices for the first and second dimension, if we append two such lattices `L1` and `L2` along the dimension `1`, the merging would be operated by connecting the output vertex `c=3` of `L1` to the input vertex `b=2` of `L2`. In this way, if we would like to start the renormalization procedure from, e.g., a `2x2` Lieb lattice (12 vertices), we can just double `L` once along both directions:
```julia
L = union(L,L,1)
L = union(L,L,2)
```
We can then build a `NdSystem` from it by simply using `NdSystem(L, H, t, lHt, J)` or any other `NdSystem` constructor.

For less trivial structures, it can be useful to first build a lattice from merging small unit cells up to some desired size and only then make the lattice periodic. For instance, a periodic honeycomb lattice could be made as follows:
```julia
using LightGraphs

D = 2                # #dimensions
a = 1; b = 2;        # vertex names
shape = (1, 1)       # shape of the unit cell
pbc = false          # Boundary conditions
Vin_1  = [a]         # Input vertices in direction d=1
Vout_1 = [b]         # Output vertices in direction d=1
Vin_2  = [a]         # Input vertices in direction d=2
Vout_2 = [b]         # Output vertices in direction d=2
unit_cell = PathGraph(2)

L = NdLattice{D}(unit_cell, shape, (Vin_1, Vin_2), (Vout_1, Vout_2), pbc)
L = union(L,L,1)
Lhex = union(L,L,2)           # One hexagon
Lhex = union(Lhex,Lhex,1)     # Three hexagons
Lhex = union(Lhex,Lhex,2)     # Nine hexagons

Lhex_pbc = pbc_from_obc(Lhex) # Build a periodic copy of it
```

## Hamiltonian and dissipators

Hamiltonian and dissipators can be built from a lattice and some input operators.
```julia
hamiltonian(L, lH, t, lHt)
```
Generate a Hamiltonian for a lattice from some local Hamiltonian `lH` and a local hopping operator `lHt` with associated nearest-neighbors coupling rate `t` as <img src="https://latex.codecogs.com/gif.latex?{\textstyle&space;\sum_i}\hat{\texttt{lH}}_i&space;&plus;&space;{\textstyle&space;\sum_{\langle&space;i;j\rangle}}&space;(t\times\hat{\texttt{lHt}}_i\otimes\hat{\texttt{lHt}}_j^\dagger&space;&plus;&space;\mathrm{H.c.})" title="{\textstyle \sum_i}\hat{\texttt{lH}}_i + {\textstyle \sum_{\langle i;j\rangle}} (t\times\hat{\texttt{lHt}}_i\otimes\hat{\texttt{lHt}}_j^\dagger + \mathrm{H.c.})" />.
```julia
hamiltonian(L, lH, (t1,t2, ... ,tN), lHt)
```
Generate a Hamiltonian for a `NdLattice` in the same way but with anisotropic nearest-neighbors coupling rates.
```julia
dissipators(L, J)
```
Generate a vector of jump operators for the lattice from a vector of local jump operators.

## Systems

Corner merging and evolving functions take systems as arguments. Systems contain a lattice, a Hamiltonian, a set of jump operators, some other operators required during merging operations and a set of observables to be transformed into the corner subspace along with the Liouvillian.

####  NdSystem
```julia
NdSystem(L, H, t, lHt, J)
NdSystem(L, H, (t1,t2, ... ,tN), lHt, J)
NdSystem(L, H, t, lHt, J, obs)
NdSystem(L, H, (t1,t2, ... ,tN), lHt, J, obs)
```
Generate a system from a `NdLattice`, a Hamiltonian, either a single nearest-neighbors coupling rate `t` (isotropic system) or a tuple with as many rates as lattice dimensions (anisotropic system), and a local hopping operator `lHt` and a vector `J` of local jump operators. If passed, `obs` contains the observables to be transformed into the corner subspace along with the Liouvillian. Either a vector of `Dict` mapping some operator names (`String`) to some global operators or a vector of `Dict` mapping some operator names (`String`) to some local operators can be passed as argument.

## Corner operations

```julia
s1Us2 = merge(s1, s2, d, ρ1, ρ2, M)
s1Us2, ρ1xρ2 = merge(s1, s2, d, ρ1, ρ2, M; return_dm=true)
```
Projection into the corner space spanned by the `M` most probable product states of the merging of two N-dimensional systems along the direction `d`. The product state in the corner basis is optionally returned.

## Steady state evaluation

```julia
ρ = steadystate_bicg(H, J, l; log=false, kwargs...)
ρ, log = steadystate_bicg(H, J, l; log=true, kwargs...)
ρ = steadystate_bicg(s, l; log=false, kwargs...)
ρ, log = steadystate_bicg(s, l; log=true, kwargs...)
```
Evaluates the steady state by solving iteratively the linear system <img src="https://latex.codecogs.com/gif.latex?\mathcal{L}\hat{\rho}&space;=&space;0" title="\mathcal{L}\hat{\rho} = 0" /> via the stabilized biconjugate gradient method with `l` `GMRES` steps. The first line of the Liouvillian is overwritten to enforce a non-trivial trace one solution, this approximation yields an error of the order of the inverse of the square of the size of the Hilbert space.

```julia
ρ = steadystate_bicg!(ρ0, H, J, l; log=false, kwargs...)
ρ, log = steadystate_bicg!(ρ0, H, J, l; log=true, kwargs...)
ρ = steadystate_bicg!(ρ0, s, l; log=false, kwargs...)
ρ, log = steadystate_bicg!(ρ0, s, l; log=true, kwargs...)
```
Same as the above but starting the iterative process from the provided density matrix.

```julia
ρ = CornerSpaceRenorm.steadystate_bicg_LtL(H, J, l; log=false, kwargs...)
ρ, log = CornerSpaceRenorm.steadystate_bicg_LtL(H, J, l; log=true, kwargs...)
ρ = CornerSpaceRenorm.steadystate_bicg_LtL(s, l; log=false, kwargs...)
ρ, log = CornerSpaceRenorm.steadystate_bicg_LtL(s, l; log=true, kwargs...)
```
Evaluates the steady state by solving iteratively the linear system <img src="https://latex.codecogs.com/gif.latex?\mathcal{L}^\dagger\mathcal{L}\vec{\rho}&plus;\vec{\mathrm{tr}}\otimes\vec{\mathrm{tr}}\cdot\vec{\rho}&space;=&space;\vec{\mathrm{tr}}" title="\mathcal{L}^\dagger\mathcal{L}\vec{\rho}+\vec{\mathrm{tr}}\otimes\vec{\mathrm{tr}}\cdot\vec{\rho} = \vec{\mathrm{tr}}" /> via the stabilized biconjugate gradient method with `l` `GMRES` steps. No approximation of the Liovillian is made in order to enfore the trace one but in practice convergence is slower and poorer.

## Plotting

Methods for `Lattice` and `AbstractSystem` are added to function `GraphPlot.gplot` of the package [GraphPlot.jl](https://github.com/JuliaGraphs/GraphPlot.jl).
```julia
gplot(L; kwargs...)
gplot(L, locs_x_in, locs_y_in; kwargs...)
gplot(s; kwargs...)
gplot(s, locs_x_in, locs_y_in; kwargs...)
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
γ = 1.   # Dissipation rate
g = 0.5γ # Local energies
V = 2γ   # Coupling rate

# 1D lattice with periodic boundary conditions. Dimension is guessed from the
# length of the shape tuple.
L = NdLattice((8,); periodic=true)

H = hamiltonian(L, g/2 * sx, V/4/2., sz)
J = dissipators(L, [sqrt(γ) * sm])

# Generate a system from a lattice, a Hamiltonian,
# a local tunnelling operator and jump operators
s = NdSystem(L, H, V/4/2., sz, J, lobs)
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
julia> Purity = 0.9740074747395815
julia> #states = 1.117951278975922
julia> <sx> = 0.2361328965173251
julia> <sy> = 0.060883828730663475
julia> <sz> = -0.9668530523327893
```
Increasing the corner dimension from `500` to `700` yields
```julia-repl
julia> Purity = 0.9747381587498526
julia> #states = 1.1149482217023885
julia> <sx> = 0.2361340393642682
julia> <sy> = 0.0623753386362691
julia> <sz> = -0.9668528544553675
```
The corner has reached convergence up to the chosen tolerance.

#### 4x4 hardcore-boson driven-dissipative lattice

The data compiled in Table 1 of Ref. <b id="f1">[[1](#f1)]</b> obtained by using an independent MATLAB implementation of the corner space renormalization procedure can be easily reproduced with the following code:

```julia
using CornerSpaceRenorm, QuantumOptics

# Local hard-core boson basis
b = FockBasis(1)

# Build some operators in the local basis
lN = number(b)
a = destroy(b)
# Build a set of local observables in the local basis
lobs = Dict("a"=>a)

# Same parameters as for Table 1 of Ref. [1]
γ  = 1. # Dissipation rate
F  = 2. # Driving strength
J  = 1. # tunneling rate
Δω = 5. # Detuning
z  = 4. # Coordination number

M = 200 # Corner Space Dimension

# 2x4 Bose-Hubbard lattice with periodic boundary conditions.
L = NdLattice((2,4); periodic=true)
# Local Hamiltonian
lH = -Δω*lN + F*(a + dagger(a))
H = hamiltonian(L, lH, -J/z, a)
C = dissipators(L, [sqrt(γ) * a])

# Generate a system from a lattice, a Hamiltonian,
# a local tunnelling operator and jump operators
s = NdSystem(L, H, -J/z, a, C, lobs)
# Compute the steady state (by brute-force integration)
ρ1 = steadystate_bicg(s, 6; tol=1e-6, verbose=true)
# Merge two systems into some corner subspace spanned by M kets
s2 = merge(s,s,1,ρ1,ρ1,M) # 4x4-site lattice
# Find the steady-state density matrix by minimizing Lρ via the BiCGStab(l=4)
# iterative method.
ρ2 = steadystate_bicg(s2, 6; tol=1e-6, verbose=true)
hermitianize!(ρ2)

# Compute the density operators in the corner space
gN = [dagger(s2.observables[i]["a"]) * s2.observables[i]["a"] for i in 1:nv(s2.lattice)]
# Compute some observables
N = [real.(expect(gN[i],ρ2)) for i in 1:nv(s2.lattice)]
mean_N = sum(N) / nv(s2.lattice)
mean_coh_re = real(sum([expect(s2.observables[i]["a"],ρ2) for i in 1:nv(s2.lattice)]) / nv(s2.lattice))
mean_g2 = sum([real(expect(gN[e.src] * gN[e.dst],ρ2) / N[e.src] / N[e.dst]) for e in edges(s2.lattice)]) / ne(s2.lattice)

# Print a line of Table 1
println("M = ",M,"\tn = ",mean_N,"\tRe(<a>) = ",mean_coh_re,"\tg2<j,l> = ",mean_g2)
```

yielding

```julia-repl
julia> M = 200 n = 0.0953582079216167  Re(<a>) = 0.27677860030342605 g2<j,l> = 1.0609718397305627
```

which indeed matches the values of the 4th row of Table 1.

#### 12x12 hardcore-boson driven-dissipative lattice

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

# Start from a 3x3 Bose-Hubbard lattice with periodic boundary conditions.
# Dimension is guessed from the length of the shape tuple.
L = NdLattice((3,3); periodic=true)
# Local Hamiltonian
lH = -Δ*N + F*(a + dagger(a))
H = hamiltonian(L, lH, -t/2., a)
J = dissipators(L, [sqrt(κ) * a])

mean_pop(s,ρ) = real(sum([expect(s.observables[i]["N"],ρ) for i in 1:nv(s.lattice)]) / nv(s.lattice))

# Generate a system from a lattice, a Hamiltonian,
# a local tunnelling operator and jump operators
s = NdSystem(L, H, -t/2., a, J, lobs)
# Find the steady-state density matrix by minimizing Lρ via the BiCGStab(l=4)
# iterative method.
@time ρ = steadystate_bicg(s, 6; tol=1e-4, verbose=true)
println("\nIteration 1")
println("Purity = ",real(tr(ρ*ρ)),"\n#states = ", real(exp(entropy_vn(ρ))))
println("<N> = ",mean_pop(s,ρ))

# Merge two systems into some corner subspace spanned by 500 kets
s = merge(s,s,1,ρ,ρ,1000) # 6x3-site lattice
# Find the steady-state density matrix
@time ρ = steadystate_bicg(s, 6; tol=1e-4, verbose=true)
hermitianize!(ρ)
println("\nIteration 2")
println("Purity = ",real(tr(ρ*ρ)),"\n#states = ", real(exp(entropy_vn(ρ))))
println("<N> = ",mean_pop(s,ρ))

# Repeat the process
s = merge(s,s,2,ρ,ρ,1000) # 6x6-site lattice
@time ρ = steadystate_bicg(s, 6; tol=1e-4, verbose=true)
hermitianize!(ρ)
println("\nIteration 3")
println("Purity = ",real(tr(ρ*ρ)),"\n#states = ", real(exp(entropy_vn(ρ))))
println("<N> = ",mean_pop(s,ρ))

# Again
s = merge(s,s,1,ρ,ρ,1000) # 12x6-site lattice
@time ρ = steadystate_bicg(s, 6; tol=1e-4, verbose=true)
hermitianize!(ρ)
println("\nIteration 4")
println("Purity = ",real(tr(ρ*ρ)),"\n#states = ", real(exp(entropy_vn(ρ))))
println("<N> = ",mean_pop(s,ρ))

# And again...
@time s = merge(s,s,2,ρ,ρ,1000) # 12x12-site lattice
@time ρ = steadystate_bicg(s, 6; tol=1e-4, verbose=true)
hermitianize!(ρ)
println("\nIteration 5")
println("Purity = ",real(tr(ρ*ρ)),"\n#states = ", real(exp(entropy_vn(ρ))))
println("<N> = ",mean_pop(s,ρ))
```
```julia-repl
julia>
# Steady-state determination time
33.391954 seconds (8.91 M allocations: 3.878 GiB, 2.59% gc time)
Iteration 1
Purity = 0.9976260308338634
#states = 1.0116079830412088
<N> = 0.008365652074378548

# Steady-state determination time
481.051435 seconds (2.01 M allocations: 17.274 GiB, 0.49% gc time)
Iteration 2
Purity = 0.8243971487739387
#states = 2.0868584800171677
<N> = 0.022055165597897344

# Steady-state determination time
1127.749763 seconds (2.98 k allocations: 21.562 GiB, 0.25% gc time)
Iteration 3
Purity = 0.647720358624363
#states = 4.140104941776228
<N> = 0.02456904145500847

# Steady-state determination time
2758.450997 seconds (553.93 k allocations: 26.773 GiB, 0.13% gc time)
Iteration 4
Purity = 0.30336292598866943
#states = 22.89421659581906
<N> = 0.03238841220844765

# Merging time
444.331974 seconds (4.01 M allocations: 9.778 GiB, 0.56% gc time)
# Steady-state determination time
5405.893048 seconds (1.00 M allocations: 30.020 GiB, 0.08% gc time)
Iteration 4
Purity = 0.02452649425569945
#states = 246.1450221114231
<N> = 0.05323380908335849
```
In the last two steps, the entropy soared, a higher corner dimensions is thus to be chosen for convergence.

## References

<b id="f1">[1]</b> S. Finazzi, A. Le Boité, F. Storme, A. Baksic, and C. Ciuti, Corner-Space Renormalization Method for Driven-Dissipative Two-Dimensional Correlated Systems, [Phys. Rev. Lett. 115, 080604 (2015)](https://doi.org/10.1103/PhysRevLett.115.080604) (thereby results were obtained by the authors using a former independent MATLAB implementation)[↩](#a1)

<b id="f2">[2]</b> W. Casteels, R. Rota, F. Storme, and C. Ciuti, Probing photon correlations in the dark sites of geometrically frustrated cavity lattices, [Phys. Rev. A. 93, 043833 (2016)](https://doi.org/10.1103/PhysRevA.93.043833)[↩](#a2)
