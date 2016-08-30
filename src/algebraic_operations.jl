function Base.:+{T<:CoordinateVector}(A::T,B::T)
    A.coord == B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    return T(A.data .+ B.data, A.coord, A.J, A.g, A.h)
end

function Base.:+{T<:CoordinateVector,S<:CoordinateVector}(A::T,B::S)
    return A + convert(B)
end

function Base.:+{T<:Real,S<:CoordinateVector}(a::T,x::S)
    return S(x.data+a,x.coord,x.J,x.g,x.h)
end

function Base.:+{T<:CoordinateVector,S<:Real}(x::T,a::S)
    return T(x.data + a,x.coord,x.J,x.g,x.h)
end

function Base.:-{T<:CoordinateVector}(A::T,B::T)
    A.coord == B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    return T(A.data .- B.data, A.coord, A.J, A.g, A.h)
end

function Base.:-{T<:CoordinateVector,S<:CoordinateVector}(A::T,B::S)
    return A - convert(B)
end

function Base.:-{T<:Real,S<:CoordinateVector}(a::T,x::S)
    return S(a - x.data,x.coord,x.J,x.g,x.h)
end

function Base.:-{T<:CoordinateVector,S<:Real}(x::T,a::S)
    return T(x.data - a,x.coord,x.J,x.g,x.h)
end

function Base.:*{T<:Real,S<:CoordinateVector}(a::T,x::S)
    return S(a*x.data,x.coord,x.J,x.g,x.h)
end

function Base.:*{T<:CoordinateVector,S<:Real}(x::T,a::S)
    return T(a*x.data,x.coord,x.J,x.g,x.h)
end

function Base.:/{T<:Real,S<:CoordinateVector}(a::T,x::S)
    return S(a./ x.data,x.coord,x.J,x.g,x.h)
end

function Base.:/{T<:CoordinateVector,S<:Real}(x::T,a::S)
    return T(x.data/a,x.coord,x.J,x.g,x.h)
end

function dot{T<:CoordinateVector}(A::T,B::T)
    A.coord === B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    dp = sum(A.data .* (A.g*B.data))
    return dp
end

function dot{T<:CoordinateVector,S<:CoordinateVector}(A::T,B::S)
    A.coord === B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    dp = sum(A.data .* B.data)
    return dp
end

norm{T<:CoordinateVector}(x::T) = sqrt(dot(x,x))

function cross{T<:CoordinateVector}(A::T,B::T)
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

function cross{T<:CoordinateVector,S<:CoordinateVector}(A::T,B::S)
    return cross(A,convert(B))
end
