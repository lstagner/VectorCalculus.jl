using VectorCalculus
using Base.Test

@testset "VectorCalculus" begin
    include("algebraic.jl")
    include("differential.jl")
end
