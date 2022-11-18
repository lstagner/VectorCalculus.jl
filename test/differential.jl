function grav_potential_cartesian(u)
    x,y,z = u
    return -1.0/sqrt(x^2 + y^2 + z^2)
end

function grav_potential_cylindrical(u)
    r,θ,z = u
    return -1.0/sqrt(r^2 + z^2)
end

function grav_potential_spherical(u)
    r,θ,ϕ = u
    return -1.0/r
end

function grav_field(potential::Function, u::Coordinate)
    return -∇(potential, u)
end

# Scalar Field Functions
function scalar_field(u::Coordinate)
    return u[3]*u[1]^2 + sin(u[1])*u[2]*u[3]^2
end

function grad_scalar_field_cartesian(u::Coordinate)
    x,y,z = u
    return [2*x*z + cos(x)*y*z^2, sin(x)*z^2, x^2 + 2*sin(x)*y*z]
end

function laplacian_scalar_field_cartesian(u::Coordinate)
    x,y,z = u
    return (2*z - sin(x)*y*z^2) + (0.0) + (2*sin(x)*y)
end

function grad_scalar_field_cylindrical(u::Coordinate)
    r,θ,z = u
    return [2*r*z + cos(r)*θ*z^2, (sin(r)*z^2)/r, r^2 + 2*sin(r)*θ*z]
end

function laplacian_scalar_field_cylindrical(u::Coordinate)
    r,θ,z = u
    return (4*z + (θ*z^2)*(cos(r)/r - sin(r))) + 2*sin(r)*θ
end

function grad_scalar_field_spherical(u::Coordinate)
    r,θ,ϕ = u
    return [2*r*ϕ + cos(r)*θ*ϕ^2, (sin(r)*ϕ^2)/(r*sin(ϕ)), (r^2 + 2*sin(r)*θ*ϕ)/r]
end

function laplacian_scalar_field_spherical(u::Coordinate)
    r,θ,ϕ = u
    return (6*ϕ + (θ*ϕ^2)*(2*cos(r)/r - sin(r))) + (cot(ϕ) + (2*θ*sin(r)/(r^2))*(cot(ϕ)*ϕ + 1))
end

#Vector Field Functions
function vector_field(u)
    V = [u[2]*u[1]^2, u[1]*u[2]^2, u[1]*u[2]*u[3]]
    return ContravariantVector(V, u, unit_basis=true)
end

function div_vector_field_cartesian(u)
    x,y,z = u
    return 5*x*y
end

function curl_vector_field_cartesian(u)
    x,y,z = u
    return [x*z, -y*z, y^2 - x^2]
end

function vector_laplacian_cartesian(u)
    x,y,z = u
    return [2*y , 2*x, 0.0]
end

function div_vector_field_cylindrical(u)
    r,θ,z = u
    return 4*θ*r + 2*θ
end

function curl_vector_field_cylindrical(u)
    r,θ,z = u
    return [z, -θ*z, 2*θ^2 - r]
end

function vector_laplacian_cylindrical(u)
    r,θ,z = u
    return [ θ*(3 - 4/r), 2*(1 + 1/r), θ*z/r]
end

function div_vector_field_spherical(u)
    r,θ,ϕ = u
    return 4*θ*r + (2*θ)/sin(ϕ) + θ*(1 + ϕ*cot(ϕ))
end

function curl_vector_field_spherical(u)
    r,θ,ϕ = u
    return [cot(ϕ)*θ^2 - ϕ/sin(ϕ), 2*θ*ϕ, (r/sin(ϕ)) - 2*θ^2]
end

function vector_laplacian_spherical(u)
    r,θ,ϕ = u
    return [4*θ - (2*θ/r)*(1 + 2/sin(ϕ) + ϕ*cot(ϕ)),
            (2*θ^2)/r + (2-θ^2)/(r*sin(ϕ)^2) + 2*(1+ ϕ*cot(ϕ)/r)/sin(ϕ),
            (θ*ϕ/r)*(2 - 1/sin(ϕ)^2) + (θ*cot(ϕ)/r)*(1 - 4/sin(ϕ))]
end

function Φ₁(c)
    tan(prod(c))
end

function Φ₂(c)
    c[1]*sin(c[2]/c[3])^2
end

function differential_identities(c)

    @test ∇(x->Φ₁(x)*Φ₂(x), c) ≈ Φ₁(c)*∇(Φ₂,c) + ∇(Φ₁,c)*Φ₂(c)

    @test ∇ ⋅ (x->Φ₁(x)*vector_field(x), c) ≈ Φ₁(c) * ∇ ⋅ (vector_field,c) + vector_field(c) ⋅ ∇(Φ₁,c)

    @test ∇ × (x->Φ₁(x)*vector_field(x), c) ≈ Φ₁(c) * ∇ × (vector_field,c) + (∇(Φ₁,c) × vector_field(c))

end

function random_cartesian()
    CartesianCoordinate(3*randn(3))
end

function random_cylindrical()
    x = [0.0, 0.0, -3.0] + [3.0,2pi,6.0].*rand(3)
    CylindricalCoordinate(x)
end

function random_spherical()
    x = [0.0, 0.0, 0.0] + [3.0,2pi,pi].*rand(3)
    SphericalCoordinate(x)
end

const ubc = unit_basis_components

@testset "Differential Operations" begin
    x = CartesianCoordinate([2.0,-1.0,3.0])
    cyl = CylindricalCoordinate([sqrt(x[1]^2 + x[2]^2), atan(x[2],x[1]), x[3]])
    sph = SphericalCoordinate([norm(x), atan(x[2],x[1]), acos(x[3]/norm(x))])

    f = norm(grav_field(grav_potential_cartesian, x))
    @test f ≈ norm(grav_field(grav_potential_cylindrical, cyl))
    @test f ≈ norm(grav_field(grav_potential_spherical, sph))

    x = random_cartesian()
    cyl = random_cylindrical()
    sph = random_spherical()
    println("Testing Gradient...")
    @testset "Gradient" begin
        @test ubc(∇(scalar_field, x)) ≈ grad_scalar_field_cartesian(x)
        @test ubc(∇(scalar_field, cyl)) ≈ grad_scalar_field_cylindrical(cyl)
        @test ubc(∇(scalar_field, sph)) ≈ grad_scalar_field_spherical(sph)
    end

    println("Testing Laplacian...")
    @testset "Laplacian" begin
        @test Δ(scalar_field, x) ≈ laplacian_scalar_field_cartesian(x)
        @test Δ(scalar_field, cyl) ≈ laplacian_scalar_field_cylindrical(cyl)
        @test Δ(scalar_field, sph) ≈ laplacian_scalar_field_spherical(sph)
    end

    println("Testing Divergence...")
    @testset "Divergence" begin
        @test ∇ ⋅ (vector_field, x) ≈ div_vector_field_cartesian(x)
        @test ∇ ⋅ (vector_field, cyl) ≈ div_vector_field_cylindrical(cyl)
        @test ∇ ⋅ (vector_field, sph) ≈ div_vector_field_spherical(sph)
    end

    println("Testing Curl...")
    @testset "Curl" begin
        @test ubc(∇ × (vector_field, x)) ≈ curl_vector_field_cartesian(x)
        @test ubc(∇ × (vector_field, cyl)) ≈ curl_vector_field_cylindrical(cyl)
        @test ubc(∇ × (vector_field, sph)) ≈ curl_vector_field_spherical(sph)
    end

    println("Testing Vector Laplacian...")
    @testset "Vector Laplacian" begin
        @test ubc(∇²(vector_field, x)) ≈ vector_laplacian_cartesian(x)
        @test ubc(∇²(vector_field, cyl)) ≈ vector_laplacian_cylindrical(cyl)
        @test ubc(∇²(vector_field, sph)) ≈ vector_laplacian_spherical(sph)
    end

    println("Testing Differential Identities...")
    @testset "Differential Identities" begin
        @testset "Cartesian" begin
            differential_identities(x)
        end
        @testset "Cylindrical" begin
            differential_identities(cyl)
        end
        @testset "Spherical" begin
            differential_identities(sph)
        end
    end
end
