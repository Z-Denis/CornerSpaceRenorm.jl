# CornerSpaceRenorm.jl

[![Build Status](https://travis-ci.com/Z-Denis/CornerSpaceRenorm.jl.svg?token=XuYcpCDomapYmd2vHj9y&branch=master)](https://travis-ci.com/Z-Denis/CornerSpaceRenorm.jl)
[![codecov](https://codecov.io/gh/Z-Denis/CornerSpaceRenorm.jl/branch/master/graph/badge.svg?token=EwifsJO3ew)](https://codecov.io/gh/Z-Denis/CornerSpaceRenorm.jl)

Skeleton of what could be a Julia package for performing corner space renormalization.

## Installation
From the Pkg REPL (prefix `]`), type:
```julia-repl
(v1.1) pkg> add https://github.com/TheorieMPQ/CornerSpaceRenorm.jl
```

## Examples
```julia
using CornerSpaceRenorm, QuantumOptics

ω = 1.
κ = 1.
F = 1e-1 * κ
J = 1.

# Local hardcore-boson basis
lb = FockBasis(1)
# Square lattice with open boundary conditions
L = SquareLattice(2,2)
# Lattices can be combined: L ∪ L (vunion(L,L)), hunion(L,L)
gb = CompositeBasis([lb for i in 1:nv(L)])
# local operators of the tunnelling Hamiltonian
tH = (J * destroy(lb),create(lb))
# Build Hamiltonian for the lattice L
H = hamiltonian(L,ω*number(lb) + F * (create(lb) + destroy(lb)),[tH,dagger.(tH)])
# Build some jump operators for the lattice L
J = dissipators(L,[sqrt(κ) * destroy(lb)])

# Generate a system from a lattice, a Hamiltonian,
# a local tunnelling operator and jump operators
s1 = System(L,H,destroy(lb),J)
# Compute the steady state (by brute-force integration)
ρ1 = steadystate.master(s1.H,s1.J)[2][end]
# Vertically merge two systems into some corner subspace
# spanned by 100 kets
s2 = vmerge(s1,s1,ρ1,ρ1,100) # 4 x 2 lattice
ρ2 = steadystate.master(s2.H,s2.J)[2][end]
# Repeat
s3 = hmerge(s2,s2,ρ2,ρ2,100) # 4 x 4 lattice
ρ3 = steadystate.master(s3.H,s3.J)[2][end]
# Again...
s4 = vmerge(s3,s3,ρ3,ρ3,100) # 8 x 4 lattice
ρ4 = steadystate.master(s4.H,s4.J)[2][end]
# etc.
s5 = hmerge(s4,s4,ρ4,ρ4,100) # 8 x 8 lattice
ρ5 = steadystate.master(s5.H,s5.J)[2][end]
println("Purity = ",real(tr(ρ5*ρ5)),"\n#states = ", real(exp(entropy_vn(ρ5))))
```
```julia-repl
julia> Purity = 0.9908348764841518
julia> #states = 1.0494305228388006
```
Doubling the corner dimension from `100` to `200` yields
```julia-repl
julia> Purity = 0.9937208692824874
julia> #states = 1.0366056109469342
```
A corner dimension of `500` yields
```julia-repl
julia> Purity = 0.9939700852528554
julia> #states = 1.0368529902181345
```
that is, from `200` to `500` the relative variation of the entropy is `0.66%` so
that we do not increase significantly our knowledge about the steady state by
going to higher corner dimensions, i.e. the corner has reached convergence.

## More generic Lattices

More generic (![alt text](https://latex.codecogs.com/gif.latex?\mathbb{Z}^N)) lattices with an arbitrary number of directions `N` along which they
can be appended are implemented via the type `ZnLattice{N}`. Its constructor
interface `ZnSystem(lat, H, lHt, J)` builds by default an `N`-dimensional cubic
lattice, however it is possible to define more generic unit cells by directly
using the default constructor. For instance, a Lieb lattice can be built as a
two-dimensional `ZnLattice` starting from a 3-vertex path graph and specifying
custom input and output vertices for each merging direction:
```julia
g = PathGraph(3)
L = ZnLattice{2}(g,(2,2),([2],[2]),([3],[1]))
```
Here `g` is a L shaped unit cell with vertices `1` (North), `2` (center) and `3`
(East). `(2, 2)` is a `2-Tuple` indicating that we embed this unit into a `2x2`
unit cell. `([2],[2])` (resp. `([3],[1])`) declares the input (resp. output)
vertices for the first and second dimension, i.e. if we append two such lattices
`L1` and `L2` along the dimension `1`, the merging would be operated by connecting
the output vertex `3` of `L1` to the input vertex `2` of `L2`. In this way, if we
would like to start the renormalization procedure from, e.g., a `2x2` Lieb lattice
(12 vertices), we can just double `L` once along both directions:
```julia
L = union(L,L,1)
L = union(L,L,2)
```
Given the operators of the `SquareLattice` example, we can then build a `ZnSystem`
from it by simply using `ZnSystem(L,H,destroy(lb),J)`.

## Plotting

Methods for `Lattice` and `AbstractSystem` are added to function `GraphPlot.gplot`
of the package [GraphPlot.jl](https://github.com/JuliaGraphs/GraphPlot.jl).
