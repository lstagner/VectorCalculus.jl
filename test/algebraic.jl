function vector_identities(c)
    A = ContravariantVector(2*randn(3),c,unit_basis=true)
    B = ContravariantVector(2*randn(3),c,unit_basis=true)
    C = CovariantVector(2*randn(3),c,unit_basis=true)
    D = CovariantVector(2*randn(3),c,unit_basis=true)

    ϵ = eps(abs(float(one(eltype(A)))))

    A_BxC = dot(A,cross(B,C))
    @test isapprox(A_BxC, dot(cross(A,B),C), atol=500ϵ)
    @test isapprox(A_BxC, dot(B,cross(C,A)), atol=500ϵ)
    @test isapprox(A_BxC, dot(cross(B,C),A), atol=500ϵ)
    @test isapprox(A_BxC, dot(C,cross(A,B)), atol=500ϵ)
    @test isapprox(A_BxC, dot(cross(C,A),B), atol=500ϵ)

    AxBxC = cross(A,cross(B,C))
    @test isapprox(AxBxC, cross(cross(C,B),A), atol=500ϵ)
    @test isapprox(AxBxC, dot(A,C)*B - dot(A,B)*C, atol=500ϵ)

    @test isapprox(cross(A,cross(B,C)) + cross(B,cross(C,A)) + cross(C,cross(A,B)), zeros(eltype(A),length(A)), atol=2000ϵ)

    AxB_CxD = dot(cross(A,B),cross(C,D))
    @test isapprox(AxB_CxD, dot(A,C)*dot(B,D) - dot(A,D)*dot(B,C), atol=2000ϵ)

    AxBxCxD = cross(cross(A,B),cross(C,D))
    @test isapprox(AxBxCxD, (dot(cross(A,B),D)*C - dot(cross(A,B),C)*D), atol=2000ϵ)
end

@testset "Algebraic Operations" begin
    println("Testing Algebraic Identities...")
    tic()
    @testset "Algebraic Identities" begin
        @testset "Cartesian" begin
            c = CartesianCoordinate([2.0,-1.0,1.5])
            vector_identities(c);
        end
        @testset "Cylindrical" begin
            c = CylindricalCoordinate([2.0,pi/4,-1.0])
            vector_identities(c);
        end
        @testset "Spherical" begin
            c = SphericalCoordinate([2.0,pi/4,pi/3])
            vector_identities(c);
        end
    end
    toc()
end
