function dot{T<:CoordinateVector}(A::T,B::T)
    A.coord=== B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    dp = (A.data./A.h)'*A.g*(B.data./B.h)
    return dp[1]
end

function dot{T<:CoordinateVector,S<:CoordinateVector}(A::T,B::S)
    A.coord=== B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    dp = (A.data./A.h)'*(B.data./B.h)
    return dp[1]
end

norm{T<:CoordinateVector}(x::T) = sqrt(dot(x,x))

using Combinatorics: levicivita
function cross{T<:CoordinateVector}(A::T,B::T)
    A.coord === B.coord || throw(ArgumentError("Vectors must have the same coordinates"))
    length(A) == 3 || throw(DomainError("Cross product only defined for 3-vectors"))
    x = zeros(eltype(A),3)

    detJ = abs(det(A.J))

    for i=1:3
        for j=1:3
            v = detJ*A[i]*B[j]/(A.h[i]*B.h[j])
            for k=1:3
                x[k] = x[k] + levicivita([i,j,k])*v
            end
        end
    end

    return T(SVector{3}(x)./A.h,A.coord,A.J,A.g,A.h)
end

function cross{T<:CoordinateVector,S<:CoordinateVector}(A::T,B::S)
    throw(ArgumentError("Vectors must have the same basis"))
end
