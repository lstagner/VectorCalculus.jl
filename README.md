# VectorCalculus

[![Build Status](https://travis-ci.org/lstagner/VectorCalculus.jl.svg?branch=master)](https://travis-ci.org/lstagner/VectorCalculus.jl)

Automatic Vector Calculus calculations in Julia via Automatic Differentiation.

See the tests for full functionality

```julia
# Some vector field function F that takes a coordinate u and returns a vector
function VF(u)
    V = [u[2]*u[1]^2, u[1]*u[2]^2, u[1]*u[2]*u[3]]
    return ContravariantVector(V, u, unit_basis=true)
end

# A scalar function that takes in a coordinate and returns a scalar
function SF(u::Coordinate)
    return u[3]*u[1]^2 + sin(u[1])*u[2]*u[3]^2
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

    @assert ∇ ⋅ (x->Φ₁(x)*vector_field(x), c) ≈ Φ₁(c) * ∇ ⋅ (vector_field,c) + vector_field(c) ⋅ ∇(Φ₁,c)

    @assert ∇ × (x->Φ₁(x)*vector_field(x), c) ≈ Φ₁(c) * ∇ × (vector_field,c) + (∇(Φ₁,c) × vector_field(c))
end

differential_identities(sph)

```
