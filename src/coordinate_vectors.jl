function metric{T}(J::AbstractArray{T,2})
    N = size(J,2)
    g = T[dot(J[:,i],J[:,j]) for i=1:N, j=1:N]
    return g
end

function covariant_basis{T}(c::Coordinate{T})
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

abstract CoordinateVector{T} <: AbstractArray{T,1}

type ContravariantVector{T,F} <: CoordinateVector{T}
    data::Vector{T}       #Contravariant Components of the Vector
    coord::Coordinate{T,F}  #Vector Coordinate
    J::Matrix{T}
    g::Matrix{T}
    h::Vector{T}
end

function ContravariantVector{T,F}(data::AbstractVector{T}, coord::Coordinate{T,F}; unit_basis = false)
    J = covariant_basis(coord)
    g = metric(J)
    h = sqrt.(diag(g))
    if unit_basis
        data = data./h
    end
    ContravariantVector{T,F}(data, coord, J, g, h)
end

type CovariantVector{T,F} <: CoordinateVector{T}
    data::Vector{T}       #Covariant Components of the Vector
    coord::Coordinate{T,F}  #Vector Coordinate
    J::Matrix{T}
    g::Matrix{T}
    h::Vector{T}
end

function CovariantVector{T,F}(data::AbstractVector{T}, coord::Coordinate{T,F}; unit_basis = false)
    J = contravariant_basis(coord)
    g = metric(J)
    h = sqrt.(diag(g))
    if unit_basis
        data = data./h
    end
    CovariantVector{T,F}(data, coord, J, g, h)
end

# Array Interface for CoordinateVector
size{T<:CoordinateVector}(A::T) = size(A.data)
getindex{T<:CoordinateVector}(A::T, i::Int) = A.data[i]
setindex!{T}(A::CoordinateVector{T}, v::T, i::Int) = A.data[i] = v
Base.linearindexing(::CoordinateVector) = Base.LinearFast()

function convert{T<:CoordinateVector,S<:CoordinateVector}(::Type{T}, x::S)
    data = x.g*x.data
    g = inv(x.g)
    J = x.J*g
    h = sqrt.(diag(g))
    T(data, x.coord, J, g, h)
end

convert{T<:CoordinateVector}(::Type{T}, x::T) = x
convert(x::CovariantVector) = convert(ContravariantVector,x)
convert(x::ContravariantVector) = convert(CovariantVector,x)

function unit_basis_components{T<:CoordinateVector}(x::T)
    return x .* x.h
end
