module VectorCalculus

using ForwardDiff
using StaticArrays

import Base: dot, cross, norm, getindex, convert, size

include("coordinates.jl")
include("coordinate_vectors.jl")
include("algebraic_operations.jl")
include("differential_operations.jl")

export CartesianCoordinate, PolarCoordinate, CylindricalCoordinate, SphericalCoordinate
export Coordinate, CoordinateVector
export CovariantVector, ContravariantVector
export metric
export covariant_basis, covariant_metric
export contravariant_basis, contravariant_metric
end
