function vector_identities(c)
    A = ContravariantVector(2*randn(3),c,unit_basis=true)
    B = ContravariantVector(2*randn(3),c,unit_basis=true)
    C = CovariantVector(2*randn(3),c,unit_basis=true)
    D = CovariantVector(2*randn(3),c,unit_basis=true)

    ϵ = eps(abs(float(one(eltype(A)))))

    A_BxC = A ⋅ (B × C)
    @test A_BxC ≈ (A × B) ⋅ C  atol=500ϵ
    @test A_BxC ≈ B ⋅ (C × A)   atol=500ϵ
    @test A_BxC ≈ (B × C) ⋅ A   atol=500ϵ
    @test A_BxC ≈ C ⋅ (A × B)   atol=500ϵ
    @test A_BxC ≈ (C × A) ⋅ B   atol=500ϵ

    AxBxC = cross(A,cross(B,C))
    AxBxC = A × (B × C)
    @test AxBxC ≈ (C × B) × A atol=500ϵ
    @test AxBxC ≈ (A ⋅ C)*B - (A ⋅ B)*C atol=500ϵ

    @test A × (B × C) + B × (C × A) + C × (A × B) ≈ zeros(eltype(A),length(A)) atol=2000ϵ

    @test (A × B) ⋅ (C × D) ≈ (A ⋅ C)*(B ⋅ D) - (A ⋅ D)*(B ⋅ C) atol=2000ϵ

    @test (A × B) × (C × D) ≈ ((A × B) ⋅ D)*C - ((A × B) ⋅ C)*D atol=2000ϵ
end

@testset "Algebraic Operations" begin
    println("Testing Algebraic Identities...")
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
end
