module VectorCalculus

using ForwardDiff
using LinearAlgebra

import LinearAlgebra: dot, cross, norm, div
import Base: div
import Base: div, getindex, setindex!, convert, size, similar, parent, axes, IndexStyle

using Base: @propagate_inbounds

using Combinatorics: levicivita

include("coordinates.jl")
include("coordinate_vectors.jl")
include("algebraic_operations.jl")
include("differential_operations.jl")

export Cartesian, Cylindrical, Spherical
export CartesianCoordinate, CylindricalCoordinate, SphericalCoordinate
export Coordinate, CoordinateVector
export CovariantVector, ContravariantVector
export metric
export covariant_basis, covariant_metric
export contravariant_basis, contravariant_metric
export unit_basis_components
export grad, laplacian, curl, vector_laplacian

end
