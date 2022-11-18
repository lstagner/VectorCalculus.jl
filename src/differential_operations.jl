struct DelOperator{T}
    x::T
end

const ∇ = DelOperator

Base.:*(::Type{∇}, x::T) where T = DelOperator{T}(x)
Base.:*(x::T, ::Type{∇}) where T = DelOperator{T}(x)

function grad(Φ::Function, u::Coordinate)
    components = ForwardDiff.gradient(Φ, u)
    return CovariantVector(components, u)
end

∇(Φ::Function, u::Coordinate) = grad(Φ,u)

function div(A::Function, u::Coordinate)
    detJ = det(covariant_basis(u))
    function detJA(x)
        a = convert(ContravariantVector, A(x))
        det(a.J) .* a
    end
    return sum(diag(ForwardDiff.jacobian(detJA, u)))./detJ
end

⋅(::Type{∇},args::Tuple{Function,Coordinate}) = div(args...)
⋅(D::∇,args::Tuple{Function,Coordinate}) = D.x*div(args...)

function laplacian(Φ::Function, u::Coordinate)
    div( x -> grad(Φ,x), u)
end

Δ(Φ::Function, u::Coordinate) = laplacian(Φ, u)

function curl(A::Function, u::Coordinate)
    n = length(u)
    J = covariant_basis(u)
    detJ = det(J)
    function coA(x)
        convert(CovariantVector, A(x))
    end
    dA = ForwardDiff.jacobian(coA, u)

    curlA = zeros(eltype(u),n)
    for k=1:n
        for i=1:n
            for j=1:n
                dAij = dA[j,i]/detJ
                curlA[k] = curlA[k] .+ (dAij*levicivita([i,j,k]))
            end
        end
    end
    return ContravariantVector(curlA,u)
end

×(::Type{∇},args::Tuple{Function,Coordinate}) = curl(args...)
×(D::∇,args::Tuple{Function,Coordinate}) = D.x*curl(args...)

function vector_laplacian(A::Function, u::Coordinate)
    divA(x) = div(A, x)
    curlA(x) = curl(A, x)
    return grad(divA, u) - curl(curlA, u)
end

∇²(A::Function, u::Coordinate) = vector_laplacian(A, u)
