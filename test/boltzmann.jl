
using Unitful
using NQCDistributions
using Random: rand!

@testset "VelocityBoltzmann" begin
    d = VelocityBoltzmann(300u"K", [1, 200, 3000], (3, 2))
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test all(out .!== 0.0)

    dist = DynamicalDistribution(d, d, (3,2))
    rand(dist)
end
