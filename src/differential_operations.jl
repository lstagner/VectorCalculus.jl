struct DelOperator{T}
    x::T
end

const ∇ = DelOperator

Base.:*(::Type{∇}, x::T) where T = DelOperator{T}(x)
Base.:*(x::T, ::Type{∇}) where T = DelOperator{T}(x)

function grad(Φ::Function, u::Coordinate{N}) where N
    components = SVector{N}(Zygote.gradient(Φ, u)[1])
    return CovariantVector(components, u)
end

∇(Φ::Function, u::Coordinate) = grad(Φ,u)

function div(A::Function, u::Coordinate)
    detJ = det(covariant_basis(u))
    function detJA(x)
        a = convert(ContravariantVector, A(x))
        det(a.J) .* a
    end
    return sum(diag(Zygote.jacobian(detJA, u)[1]))./detJ
end

⋅(::Type{∇},args::Tuple{Function,Coordinate}) = div(args...)
⋅(D::∇,args::Tuple{Function,Coordinate}) = D.x*div(args...)

function laplacian(Φ::Function, u::Coordinate)
    div( x -> grad(Φ,x), u)
end

Δ(Φ::Function, u::Coordinate) = laplacian(Φ, u)

function curl(A::Function, u::Coordinate{N}) where N
    n = length(u)
    J = covariant_basis(u)
    detJ = det(J)
    function coA(x)
        convert(CovariantVector, A(x))
    end
    dA = Zygote.jacobian(coA, u)[1]

    curlA = zeros(eltype(u),n)
    for k=1:n
        for i=1:n
            for j=1:n
                dAij = dA[j,i]/detJ
                curlA[k] = curlA[k] .+ (dAij*levicivita([i,j,k]))
            end
        end
    end
    return ContravariantVector(SVector{N}(curlA),u)
end

×(::Type{∇},args::Tuple{Function,Coordinate}) = curl(args...)
×(D::∇,args::Tuple{Function,Coordinate}) = D.x*curl(args...)

function vector_laplacian(A::Function, u::Coordinate)
    divA(x) = div(A, x)
    curlA(x) = curl(A, x)
    return grad(divA, u) - curl(curlA, u)
end

∇²(A::Function, u::Coordinate) = vector_laplacian(A, u)
