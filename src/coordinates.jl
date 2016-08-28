type Coordinate{T,F<:Function} <: AbstractArray{T,1}
    u::Vector{T} #Coordinates
    R::F         #Transformation (x,y,z,...) = R(u₁,u₂,u₃,...)
end

# Array Interface for Coordinates
size(A::Coordinate) = size(A.u)
getindex(A::Coordinate, i::Int) = A.u[i]
setindex!{T,F}(A::Coordinate{T,F}, v::T, i::Int) = A.u[i] = v
Base.linearindexing(::Coordinate) = Base.LinearFast()
similar{T,F}(A::Coordinate{T,F}) = Coordinate(Array{T}(size(A)), A.R)
similar{T}(A::Coordinate, ::Type{T}) = Coordinate(Array{T}(size(A)), A.R)
function similar{T,N}(A::Coordinate, ::Type{T}, dims::Tuple{Vararg{Int64,N}})
    return Coordinate(Array{T}(dims), A.R)
end

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
