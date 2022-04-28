using Test
using NQCDistributions
using Unitful, UnitfulAtomic

@testset "PureState" begin
    d = PureState(1)
    NQCDistributions.density_matrix(d, 2)
end

@testset "MixedState" begin
    d = MixedState([0.3, 0.4])
    NQCDistributions.density_matrix(d)
end

@testset "FermiDiracState" begin
    d = FermiDiracState(0.0, 300u"K")
    NQCDistributions.density_matrix(d, austrip.((-1:0.01:1) .* u"eV"))
end
