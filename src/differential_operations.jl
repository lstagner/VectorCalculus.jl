function grad(Φ::Function, u::Coordinate)
    components = ForwardDiff.gradient(Φ, u, usecache=false)
    return CovariantVector(components, u)
end

function div(A::Function, u::Coordinate)
    sqrtg = sqrt(det(covariant_metric(u)))
    function sqrtgA(x)
        a = convert(ContravariantVector, A(x))
        sqrt(det(a.g)) .* a
    end
    return sum(diag(ForwardDiff.jacobian(sqrtgA, u, usecache=false)))./sqrtg
end

function laplacian(Φ::Function, u::Coordinate)
    div( x -> grad(Φ,x), u)
end

#function curl(A::Function, u::Coordinate)
# This is going to be calculate the vector jacobian
#end
