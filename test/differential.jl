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
    return -grad(potential, u)
end

x = CartesianCoordinate([2.0,-1.0,3.0])
cyl = CylindricalCoordinate([sqrt(x[1]^2 + x[2]^2), atan2(x[2],x[1]), x[3]])
sph = SphericalCoordinate([norm(x), atan2(x[2],x[1]), acos(x[3]/norm(x))])

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

const ubc = unit_basis_components

@testset "Differential Operations" begin
    f = norm(grav_field(grav_potential_cartesian, x))
    @test f ≈ norm(grav_field(grav_potential_cylindrical, cyl))
    @test f ≈ norm(grav_field(grav_potential_spherical, sph))

    @test ubc(grad(scalar_field, x)) ≈ grad_scalar_field_cartesian(x)
    @test ubc(grad(scalar_field, cyl)) ≈ grad_scalar_field_cylindrical(cyl)
    @test ubc(grad(scalar_field, sph)) ≈ grad_scalar_field_spherical(sph)

    @test laplacian(scalar_field, x) ≈ laplacian_scalar_field_cartesian(x)
    @test laplacian(scalar_field, cyl) ≈ laplacian_scalar_field_cylindrical(cyl)
    @test laplacian(scalar_field, sph) ≈ laplacian_scalar_field_spherical(sph)
end

