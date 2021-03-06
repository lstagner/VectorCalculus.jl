struct Coordinate{T,F<:Function} <: AbstractArray{T,1}
    u::Vector{T} #Coordinates
    R::F         #Transformation (x,y,z,...) = R(u₁,u₂,u₃,...)
end

# Array Interface for Coordinates
parent(A::T) where {T<:Coordinate} = A.u
size(A::T) where {T<:Coordinate} = size(A.u)
axes(A::T) where {T<:Coordinate} = axes(A.u)
parenttype(::Type{Coordinate{T,F}}) where {T,F} = Vector{T}
IndexStyle(::Type{T}) where {T<:Coordinate} = IndexStyle(parenttype(T))

@propagate_inbounds getindex(A::T, i::Int) where {T<:Coordinate} = A.u[i]
@propagate_inbounds setindex!(A::Coordinate{T}, v::T, i::Int) where {T} = A.u[i] = v

similar(A::Coordinate{T,F}) where {T,F} = Coordinate(Array{T}(undef, size(A)), A.R)
similar(A::Coordinate, ::Type{T}) where {T} = Coordinate(Array{T}(undef, size(A)), A.R)
function similar(A::Coordinate, ::Type{T}, dims::Tuple{Vararg{Int64,N}}) where {T,N}
    return Coordinate(Array{T}(undef,dims), A.R)
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
