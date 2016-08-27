function metric{T}(J::AbstractArray{T,2})
    N = size(J,2)
    g = T[dot(J[:,i],J[:,j]) for i=1:N, j=1:N]
    return g
end

function covariant_basis(c::Coordinate)
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

abstract CoordinateVector{D,T} <: AbstractArray{T,1}

type ContravariantVector{D,T,F} <: CoordinateVector{D,T}
    data::SVector{D,T}       #Contravariant Components of the Vector
    coord::Coordinate{D,T,F}  #Vector Coordinate
    J::SMatrix{D,D,T}
    g::SMatrix{D,D,T}
    h::SVector{D,T}
end

function ContravariantVector{D,T,F}(data::AbstractVector{T}, coord::Coordinate{D,T,F}; normed=true)
    J = covariant_basis(coord)
    g = metric(J)
    h = sqrt.(diag(g))
    if normed
        c = data
    else
        c = data.*h
    end
    ContravariantVector{D,T,F}(c,coord,SMatrix{D,D,T}(J),SMatrix{D,D,T}(g),SVector{D,T}(h))
end

type CovariantVector{D,T,F} <: CoordinateVector{D,T}
    data::SVector{D,T}       #Covariant Components of the Vector
    coord::Coordinate{D,T,F}  #Vector Coordinate
    J::SMatrix{D,D,T}
    g::SMatrix{D,D,T}
    h::SVector{D,T}
end

function CovariantVector{D,T,F}(data::AbstractVector{T}, coord::Coordinate{D,T,F}; normed=true)
    J = contravariant_basis(coord)
    g = metric(J)
    h = sqrt(diag(g))
    if normed
        c = data
    else
        c = data.*h
    end
    CovariantVector{D,T,F}(c,coord,SMatrix{D,D,T}(J),SMatrix{D,D,T}(g),SVector{D,T}(h))
end

# Array Interface for CoordinateVector
size{T<:CoordinateVector}(A::T) = size(A.data)
getindex{T<:CoordinateVector}(A::T, i::Int) = A.data[i]

function convert{T<:CoordinateVector,S<:CoordinateVector}(::Type{T},x::S)
    c = x.g*(x.data ./ x.h)
    g = inv(x.g)
    J = x.J*g
    h = SVector(sqrt.(diag(g)))
    T(c.*h,x.coord,J,g,h)
end

convert{T<:CoordinateVector}(::Type{T},x::T) = x
convert(x::CovariantVector) = convert(ContravariantVector,x)
convert(x::ContravariantVector) = convert(CovariantVector,x)
