
using Unitful
using NQCDistributions
using Random: rand!
using NQCDynamics

@testset "VelocityBoltzmann" begin
    d = VelocityBoltzmann(300u"K", [1, 200, 3000], (3, 2))
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test all(out .!== 0.0)

    @test_throws DimensionMismatch DynamicalDistribution(d, d, (3,3))
    dist = DynamicalDistribution(d, d, (3,2))
    rand(dist)
    @test_throws DimensionMismatch DynamicalDistribution(d, d, (3,3,4))
    dist = DynamicalDistribution(d, d, (3,2,4))
    rand(dist)

    simulation=Simulation(Atoms([:N, :H, :H, :H]), Free(3), 300u"K")
    @test size(VelocityBoltzmann(300u"K", simulation))==(3,4)
    rp_simulation=RingPolymerSimulation(Atoms([:N, :H, :H, :H]), Free(3), 300u"K", 4)
end
