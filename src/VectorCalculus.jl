module VectorCalculus

import Base: dot, cross, norm
import Base: getindex, setindex!, convert, size, similar

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
export unit_basis_components
export grad, div, laplacian,curl

using ForwardDiff

end
