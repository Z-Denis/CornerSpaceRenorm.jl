# CornerSpaceRenorm

[![Build Status](https://travis-ci.com/Z-Denis/CornerSpaceRenorm.jl.svg?token=XuYcpCDomapYmd2vHj9y&branch=master)](https://travis-ci.com/Z-Denis/CornerSpaceRenorm.jl)
[![codecov](https://codecov.io/gh/Z-Denis/CornerSpaceRenorm.jl/branch/master/graph/badge.svg?token=EwifsJO3ew)](https://codecov.io/gh/Z-Denis/CornerSpaceRenorm.jl)

Skeleton of what could be a Julia package for performing corner space renormalisation. It currently requires a lot of RAM for large lattices due to all external tunneling terms being stored (which number grows quadratically with the size of the edge).

## Examples
```julia
using CornerSpaceRenorm, QuantumOptics

# Local hardcore-boson basis
lb = FockBasis(1)
# Square lattice with open boundary conditions
L = SquareLattice(2,2)
# Lattices can be combined: L ∪ L (vunion(L,L)), hunion(L,L)
gb = CompositeBasis([lb for i in 1:nv(L)])
# local operators of the tunnelling Hamiltonian
tH = (destroy(lb),create(lb))
# Build Hamiltonian for the lattice L
H = hamiltonian(L,number(lb) + 1e-1 * (create(lb) + destroy(lb)),[tH,dagger.(tH)])
# Build some jump operators for the lattice L
J = dissipators(L,[destroy(lb),1e-2 * create(lb)])

# Generate a system from a lattice, a Hamiltonian,
# a local tunnelling operator and jump operators
s1 = System(L,H,destroy(lb),J)
# Compute the steady state
ρ1 = steadystate.master(s1.H,s1.J)[2][end]
# Vertically merge two systems
s2 = vmerge(s1,s1) # 4 x 2 lattice
# Compute a corner subspace spanned by 50 kets
cspace, pairs = corner_subspace(ρ1,ρ1,50)
# cornerize the doubled system
s2c = cornerize(s2,cspace)
# Repeat
ρ2 = steadystate.master(s2c.H,s2c.J)[2][end]
s3 = hmerge(s2c,s2c) # 4 x 4 lattice
cspace, pairs = corner_subspace(ρ2,ρ2,50)
s3c = cornerize(s3,cspace)
# Again...
ρ3 = steadystate.master(s3c.H,s3c.J)[2][end]
s4 = vmerge(s3c,s3c) # 8 x 4 lattice
cspace, pairs = corner_subspace(ρ3,ρ3,50)
s4c = cornerize(s4,cspace)
# etc.
ρ4 = steadystate.master(s4c.H,s4c.J)[2][end]
s5 = hmerge(s4c,s4c) # 8 x 8 lattice
cspace, pairs = corner_subspace(ρ4,ρ4,50)
s5c = cornerize(s5,cspace)
ρ5 = steadystate.master(s5c.H,s5c.J)[2][end]
println("Purity = ",real(tr(ρ5*ρ5)),"\n#states = ", real(exp(entropy_vn(ρ5))))
```
```julia-repl
julia> Purity = 0.9784337030301732
julia> #states = 1.1039352555204354
```
Doubling the dimension of the corner from `50` to `100` yields:
```julia-repl
julia> Purity = 0.9731031343429087
julia> #states = 1.1333972699083248
```
