using NQCDistributions
using Test, SafeTestsets

using Distributions: Normal
using Random: rand!
using RingPolymerArrays: RingPolymerArray, eachbead
using Unitful

@testset "UnivariateFill" begin
    d = UnivariateFill(Normal(), (3,2))
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    @test d[1] isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test all(out .!== 0.0)
end

@testset "UnivariateArray" begin
    d = UnivariateArray([Normal() for i=1:3, j=1:2])
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    @test d[1] isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test all(out .!== 0.0)

    d = UnivariateArray([Normal() for i=1:3, j=1:2, k=1:4])
    @test eltype(d) == Array{Float64,3}
    @test rand(d) isa AbstractArray
    @test d[1] isa AbstractArray
    out = zeros(3,2,4)
    rand!(out, d)
    @test all(out .!== 0.0)
end

@testset "FixedArray" begin
    d = FixedArray([1 2; 3 4; 5 6])
    @test eltype(d) == Matrix{Int}
    @test rand(d) isa AbstractMatrix
    @test d[1] isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test out == [1.0 2.0; 3.0 4.0; 5.0 6.0]

    a = rand(3,2,3)
    d = FixedArray(a)
    @test eltype(d) == Array{Float64,3}
    @test rand(d) isa AbstractArray{Float64,3}
    @test d[1] isa AbstractArray{Float64,3}
    out = zeros(3,2,3)
    rand!(out, d)
    @test out == a
end

@testset "FixedValue" begin
    d = FixedValue(1.0, (3,2))
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    @test d[1] isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test out == ones(3,2)
end

@testset "ConfigurationVector" begin
    d = ConfigurationVector([ones(3,2) for _=1:3])
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    @test d[1] isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test out == ones(3,2)

    d = ConfigurationVector([ones(3,2,3) for _=1:3])
    @test eltype(d) == Array{Float64,3}
    @test rand(d) isa AbstractArray{Float64,3}
    out = zeros(3,2,3)
    rand!(out, d)
    @test out == ones(3,2,3)
end

@testset "RingPolymerWrapper" begin
    nbeads = 4
    classical = [1]

    @testset "UnivariateFill" begin
        dsingle = UnivariateFill(Normal(), (3,2))
        d = RingPolymerWrapper(dsingle, nbeads; classical)
        @test eltype(d) <: RingPolymerArray
        @test d[1] isa RingPolymerArray
        @test rand(d) isa RingPolymerArray
    end

    @testset "UnivariateArray" begin
        dsingle = UnivariateArray([Normal() for i=1:3, j=1:2])
        d = RingPolymerWrapper(dsingle, nbeads; classical)
        @test eltype(d) <: RingPolymerArray
        @test d[1] isa RingPolymerArray
        @test rand(d) isa RingPolymerArray
    end

    @testset "FixedArray" begin
        dsingle = UnivariateArray([Normal() for i=1:3, j=1:2])
        d = RingPolymerWrapper(dsingle, nbeads; classical)
        @test eltype(d) <: RingPolymerArray
        @test d[1] isa RingPolymerArray
        @test rand(d) isa RingPolymerArray
    end

    @testset "FixedValue" begin
        dsingle = FixedValue(1.0, (3,2))
        d = RingPolymerWrapper(dsingle, nbeads; classical)
        @test eltype(d) <: RingPolymerArray
        @test d[1] isa RingPolymerArray
        @test rand(d) isa RingPolymerArray
        @test rand(d) == ones(3,2,4)
    end

    @testset "ConfigurationVector" begin
        dsingle = ConfigurationVector([ones(3,2) for _=1:3])
        d = RingPolymerWrapper(dsingle, nbeads; classical)
        @test eltype(d) <: RingPolymerArray
        @test d[1] isa RingPolymerArray
        @test rand(d) isa RingPolymerArray
        @test rand(d) == ones(3,2,4)
    end
end

@testset "VelocityBoltzmann" begin
    d = VelocityBoltzmann(300u"K", [1, 200, 3000], (3, 2))
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    @test d[1] isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test all(out .!== 0.0)
end

@safetestset "HarmonicWigner" begin include("harmonic_wigner.jl") end
@safetestset "PositionHarmonicRingPolymer" begin include("harmonic_ring_polymer.jl") end

@safetestset "ElectronicDistributions" begin include("electronic.jl") end
