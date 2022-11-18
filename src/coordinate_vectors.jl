function metric(J::AbstractArray{T,2}) where {T}
    g = J'*J
    return g
end

function covariant_basis(c::Coordinate{N,T}) where {N,T}
    return SMatrix{N,N}(Zygote.jacobian(c.R,c.u)[1])
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

abstract type CoordinateVector{N,T} <: AbstractArray{T,1} end

struct ContravariantVector{N,T} <: CoordinateVector{N,T}
    data::SVector{N,T}       #Contravariant Components of the Vector
    coord::Coordinate{N,T}  #Vector Coordinate
    J::SMatrix{N,N,T}
    g::SMatrix{N,N,T}
    h::SVector{N,T}
end

function ContravariantVector(d::AbstractVector{T}, coord::Coordinate{N,T}; unit_basis = false) where {N,T}
    data = SVector{N}(d)
    J = covariant_basis(coord)
    g = metric(J)
    h = sqrt.(diag(g))
    if unit_basis
        data = data./h
    end
    ContravariantVector(data, coord, J, g, h)
end

struct CovariantVector{N,T} <: CoordinateVector{N,T}
    data::SVector{N,T}       #Covariant Components of the Vector
    coord::Coordinate{N,T}  #Vector Coordinate
    J::SMatrix{N,N,T}
    g::SMatrix{N,N,T}
    h::SVector{N,T}
end

function CovariantVector(d::AbstractVector{T}, coord::Coordinate{N,T}; unit_basis = false) where {N,T}
    data = SVector{N}(d)
    J = contravariant_basis(coord)
    g = metric(J)
    h = sqrt.(diag(g))
    if unit_basis
        data = data./h
    end
    CovariantVector(data, coord, J, g, h)
end

using Zygote: @adjoint
@adjoint CovariantVector(d, u; kwargs...) = CovariantVector(d,u; kwargs...), c⁻ -> (c⁻.data,c⁻.coord)
@adjoint ContravariantVector(d,u; kwargs...) = ContravariantVector(d,u; kwargs...), c⁻ -> (c⁻.data,c⁻.coord)

# Array Interface for CoordinateVector
parent(A::T) where {T<:CoordinateVector} = A.data
size(A::T) where {T<:CoordinateVector} = size(A.data)
axes(A::T) where {T<:CoordinateVector} = axes(A.data)
parenttype(::Type{S}) where S <: CoordinateVector{N,T} where {N,T} = Vector{T}
IndexStyle(::Type{T}) where {T<:CoordinateVector} = IndexStyle(parenttype(T))

@propagate_inbounds getindex(A::T, i::Int) where {T<:CoordinateVector} = A.data[i]
#@propagate_inbounds setindex!(A::CoordinateVector{T}, v::T, i::Int) where {T} = A.data[i] = v

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
