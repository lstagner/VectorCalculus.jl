using Test
using Random
using VectorCalculus

Random.seed!(1234)

@testset "VectorCalculus" begin
    using LinearAlgebra
    include("algebraic.jl")
    include("differential.jl")
end
