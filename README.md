# CornerSpaceRenorm

[![Build Status](https://travis-ci.com/Z-Denis/CornerSpaceRenorm.jl.svg?token=XuYcpCDomapYmd2vHj9y&branch=master)](https://travis-ci.com/Z-Denis/CornerSpaceRenorm.jl)
[![codecov](https://codecov.io/gh/Z-Denis/CornerSpaceRenorm.jl/branch/master/graph/badge.svg?token=EwifsJO3ew)](https://codecov.io/gh/Z-Denis/CornerSpaceRenorm.jl)

Skeleton of what could be a Julia package for performing corner space renormalisation.

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
s5 = hmerge(s4c,s4c,ρ4,ρ4,100) # 8 x 8 lattice
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
