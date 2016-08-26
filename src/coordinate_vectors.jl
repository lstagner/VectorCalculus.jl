abstract CoordinateVector{D,T} <: AbstractArray{T,1}

type ContravariantVector{D,T,F} <: CoordinateVector{D,T}
    data::SVector{D,T}       #Contravariant Components of the Vector
    coord::Coordinate{D,T,F}  #Vector Coordinate
    J::SMatrix{D,D,T}
    g::SMatrix{D,D,T}
    h::SVector{D,T}
end

function ContravariantVector{D,T,F}(data::AbstractVector{T}, coord::Coordinate{D,T,F})
    J = ForwardDiff.jacobian(coord.R,coord.u)
    g = T[dot(J[:,i],J[:,j]) for i=1:D, j=1:D]
    h = sqrt(diag(g))
    c = data
    ContravariantVector{D,T,F}(data,coord,SMatrix{D,D,T}(J),SMatrix{D,D,T}(g),SVector{D,T}(h))
end

type CovariantVector{D,T,F} <: CoordinateVector{D,T}
    data::SVector{D,T}       #Covariant Components of the Vector
    coord::Coordinate{D,T,F}  #Vector Coordinate
    J::SMatrix{D,D,T}
    g::SMatrix{D,D,T}
    h::SVector{D,T}
end

function CovariantVector{D,T,F}(data::AbstractVector{T}, coord::Coordinate{D,T,F})
    J = ForwardDiff.jacobian(coord.R,coord.u)
    g = T[dot(J[:,i],J[:,j]) for i=1:D, j=1:D]
    g = inv(g)
    J = g'*J
    h = sqrt(diag(g))
    c = data
    CovariantVector{D,T,F}(data,coord,SMatrix{D,D,T}(J),SMatrix{D,D,T}(g),SVector{D,T}(h))
end

# Array Interface for CoordinateVector
size{T<:CoordinateVector}(A::T) = size(A.data)
getindex{T<:CoordinateVector}(A::T, i::Int) = A.data[i]

function convert{T<:CoordinateVector,S<:CoordinateVector}(::Type{T},x::S)
    c = x.g'*(x.data ./ x.h)
    g = inv(x.g)
    J = g'*x.J
    h = SVector(sqrt.(diag(g)))
    T(c.*h,x.coord,J,g,h)
end

convert{T<:CoordinateVector}(::Type{T},x::T) = x
convert(x::CovariantVector) = convert(ContravariantVector,x)
convert(x::ContravariantVector) = convert(CovariantVector,x)
