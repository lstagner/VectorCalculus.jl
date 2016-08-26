immutable Coordinate{D,T,F<:Function} <: AbstractArray{T,1}
    u::SVector{D,T} #Coordinates
    R::F            #Transformation (x,y,z,...) = R(u₁,u₂,u₃,...)
end

function Coordinate{T<:AbstractVector}(u::T,R::Function)
    Coordinate(SVector{length(u),eltype(u)}(u),R)
end

# Array Interface for Coordinates
size(A::Coordinate) = size(A.u)
getindex(A::Coordinate, i::Int) = A.u[i]

# Standard Transformations
function Cartesian(u)
    u
end

function CartesianCoordinate(u)
    Coordinate(u,Cartesian)
end

function Cylindrical(u)
    # u -> (r,θ,z) r ∈ [0,∞), θ ∈ [0,2π), z ∈ (-∞,∞)
    r = u[1]
    θ = u[2]
    z = u[3]
    [r*cos(θ), r*sin(θ), z]
end

function CylindricalCoordinate(u)
    Coordinate(u,Cylindrical)
end

function Spherical(u)
    # u -> (r,θ,ϕ) r ∈ [0,∞), θ ∈ [0,2π), ϕ ∈ [0,π]
    r = u[1]
    θ = u[2]
    ϕ = u[3]
    [r*cos(θ)*sin(ϕ), r*sin(θ)*sin(ϕ), r*cos(ϕ)]
end

function SphericalCoordinate(u)
    Coordinate(u,Spherical)
end
