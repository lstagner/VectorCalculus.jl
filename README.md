# VectorCalculus

[![Build Status](https://travis-ci.org/lstagner/VectorCalculus.jl.svg?branch=master)](https://travis-ci.org/lstagner/VectorCalculus.jl)

Automatic Vector Calculus calculations in Julia via Automatic Differentiation.

See the tests for full functionality

```julia
# In VectorCalculus.jl a coordinate is defined by the coordinate values plus a transform function that takes the different coordinate values and transform them into cartesian.

u = rand(3)
c = Coordinate(u, x -> [x[1]*cos(x[2]), x[1]*sin(x[2]), x[3]]) # defines a cylindrical coordinate. 
c = CylindricalCoordinate(u) # is provided for convienience

# Coordinate Vectors are vectors in the coordinate system at a given point. There are two types of Coordinate Vectors: Covariant and Contravariant.
# Each coordinate vector carries around its: data, coordinate, metric, Jacobian, and basis scale factors.
# Conversion between Co and Contravariant Vectors are handled internally.

# Some vector field function F that takes a coordinate u and returns a vector
function VF(u)
    V = [u[2]*u[1]^2, u[1]*u[2]^2, u[1]*u[2]*u[3]]
    return ContravariantVector(V, u, unit_basis=true)
end

# Two scaler functions that take a Coordinate c and returns a scaler
function Φ₁(c)
    tan(prod(c))
end

function Φ₂(c)
    c[1]*sin(c[2]/c[3])^2
end

x = CartesianCoordinate([2.0,-1.0,3.0])
cyl = CylindricalCoordinate([sqrt(x[1]^2 + x[2]^2), atan(x[2],x[1]), x[3]])
sph = SphericalCoordinate([norm(x), atan(x[2],x[1]), acos(x[3]/norm(x))])

function differential_identities(c)
    @assert ∇(x->Φ₁(x)*Φ₂(x), c) ≈ Φ₁(c)*∇(Φ₂,c) + ∇(Φ₁,c)*Φ₂(c)

    @assert ∇ ⋅ (x->Φ₁(x)*VF(x), c) ≈ Φ₁(c) * ∇ ⋅ (VF,c) + VF(c) ⋅ ∇(Φ₁,c)

    @assert ∇ × (x->Φ₁(x)*VF(x), c) ≈ Φ₁(c) * ∇ × (VF,c) + (∇(Φ₁,c) × VF(c))
end

differential_identities(sph)

```
