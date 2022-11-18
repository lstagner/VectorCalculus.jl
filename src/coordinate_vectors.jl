function metric(J::AbstractArray{T,2}) where {T}
    g = J'*J
    return g
end

function covariant_basis(c::Coordinate{T}) where {T}
    return ForwardDiff.jacobian(c.R,c.u)
end

function covariant_metric(c::Coordinate)
    J = covariant_basis(c)
    return metric(J)
end

function contravariant_basis(c::Coordinate)
    Jco = covariant_basis(c)
    gco = metric(Jco)
    g = inv(gco)
    J = Jco*g
    return J
end

function contravariant_metric(c::Coordinate)
    J = contravariant_basis(c)
    g = metric(J)
end

abstract type CoordinateVector{T} <: AbstractArray{T,1} end

struct ContravariantVector{T,F} <: CoordinateVector{T}
    data::Vector{T}       #Contravariant Components of the Vector
    coord::Coordinate{T,F}  #Vector Coordinate
    J::Matrix{T}
    g::Matrix{T}
    h::Vector{T}
end

function ContravariantVector(data::AbstractVector{T}, coord::Coordinate{T,F}; unit_basis = false) where {T,F}
    J = covariant_basis(coord)
    g = metric(J)
    h = sqrt.(diag(g))
    if unit_basis
        data = data./h
    end
    ContravariantVector{T,F}(data, coord, J, g, h)
end

struct CovariantVector{T,F} <: CoordinateVector{T}
    data::Vector{T}       #Covariant Components of the Vector
    coord::Coordinate{T,F}  #Vector Coordinate
    J::Matrix{T}
    g::Matrix{T}
    h::Vector{T}
end

function CovariantVector(data::AbstractVector{T}, coord::Coordinate{T,F}; unit_basis = false) where {T,F}
    J = contravariant_basis(coord)
    g = metric(J)
    h = sqrt.(diag(g))
    if unit_basis
        data = data./h
    end
    CovariantVector{T,F}(data, coord, J, g, h)
end

# Array Interface for CoordinateVector
parent(A::T) where {T<:CoordinateVector} = A.data
size(A::T) where {T<:CoordinateVector} = size(A.data)
axes(A::T) where {T<:CoordinateVector} = axes(A.data)
parenttype(::Type{S}) where S <: CoordinateVector{T} where T = Vector{T}
IndexStyle(::Type{T}) where {T<:CoordinateVector} = IndexStyle(parenttype(T))

@propagate_inbounds getindex(A::T, i::Int) where {T<:CoordinateVector} = A.data[i]
@propagate_inbounds setindex!(A::CoordinateVector{T}, v::T, i::Int) where {T} = A.data[i] = v

function convert(::Type{T}, x::S) where {T<:CoordinateVector,S<:CoordinateVector}
    data = x.g*x.data
    g = inv(x.g)
    J = x.J*g
    h = sqrt.(diag(g))
    T(data, x.coord, J, g, h)
end

convert(::Type{T}, x::T) where {T<:CoordinateVector} = x
convert(x::CovariantVector) = convert(ContravariantVector,x)
convert(x::ContravariantVector) = convert(CovariantVector,x)

function unit_basis_components(x::T) where {T<:CoordinateVector}
    return x .* x.h
end
