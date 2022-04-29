using NQCDistributions
using Distributions: Normal

@testset "DynamicalDistribution" begin
    d = DynamicalDistribution(Normal(), 2.0, (3,2))
    @test d[1] == d[1]
    @test rand(d) != rand(d)
    @test length(collect(d)) == 1
end
