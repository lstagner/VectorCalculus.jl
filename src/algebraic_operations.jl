function Base.:+(A::T,B::T) where T <: CoordinateVector
    A.coord == B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    return T(A.data .+ B.data, A.coord, A.J, A.g, A.h)
end

function Base.:+(A::T,B::S) where {T<:CoordinateVector, S<:CoordinateVector}
    return A + convert(B)
end

function Base.:+(a::T,x::S) where {T<:Real, S<:CoordinateVector}
    return S(x.data+a,x.coord,x.J,x.g,x.h)
end

function Base.:+(x::T,a::S) where {T<:CoordinateVector, S<:Real}
    return T(x.data + a,x.coord,x.J,x.g,x.h)
end

function Base.:-(A::T,B::T) where {T<:CoordinateVector}
    A.coord == B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    return T(A.data .- B.data, A.coord, A.J, A.g, A.h)
end

function Base.:-(A::T,B::S) where {T<:CoordinateVector,S<:CoordinateVector}
    return A - convert(B)
end

function Base.:-(a::T,x::S) where {T<:Real,S<:CoordinateVector}
    return S(a - x.data,x.coord,x.J,x.g,x.h)
end

function Base.:-(x::T,a::S) where {T<:CoordinateVector,S<:Real}
    return T(x.data - a,x.coord,x.J,x.g,x.h)
end

function Base.:*(a::T,x::S) where {T<:Real,S<:CoordinateVector}
    return S(a*x.data,x.coord,x.J,x.g,x.h)
end

function Base.:*(x::T,a::S) where {T<:CoordinateVector,S<:Real}
    return T(a*x.data,x.coord,x.J,x.g,x.h)
end

function Base.:/(a::T,x::S) where {T<:Real,S<:CoordinateVector}
    return S(a./ x.data,x.coord,x.J,x.g,x.h)
end

function Base.:/(x::T,a::S) where {T<:CoordinateVector,S<:Real}
    return T(x.data/a,x.coord,x.J,x.g,x.h)
end

function dot(A::T,B::T) where {T<:CoordinateVector}
    A.coord === B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    dp = sum(A.data .* (A.g*B.data))
    return dp
end

function dot(A::T,B::S) where {T<:CoordinateVector,S<:CoordinateVector}
    A.coord === B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    dp = sum(A.data .* B.data)
    return dp
end

function ⋅(A,B)
    dot(A,B)
end

norm(x::T) where {T<:CoordinateVector} = sqrt(dot(x,x))

function cross(A::T,B::T) where {T<:CoordinateVector}
    A.coord === B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    length(A) == 3 || throw(DomainError("Cross product only defined for 3-vectors"))
    AxB = zeros(eltype(A),3)
    detJ = det(A.J)

    for i=1:3
        for j=1:3
            v = detJ*A[i]*B[j]
            for k=1:3
                AxB[k] = AxB[k] + levicivita([i,j,k])*v
            end
        end
    end

    AxB = inv(A.g)*AxB
    return T(AxB,A.coord,A.J,A.g,A.h)
end

function cross(A::T,B::S) where {T<:CoordinateVector,S<:CoordinateVector}
    return cross(A,convert(B))
end

function ×(A, B)
    return cross(A,B)
end
