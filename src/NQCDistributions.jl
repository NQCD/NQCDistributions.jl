module NQCDistributions

using Distributions: Sampleable, Univariate
using Random: Random, AbstractRNG, SamplerTrivial
using RingPolymerArrays: RingPolymerArray, eachbead

export UnivariateFill
export UnivariateArray
export FixedArray
export FixedValue
export ConfigurationVector
export RingPolymerWrapper

"""
    UnivariateFill{S<:Sampleable{Univariate}}

Fill all degrees of freedom from single Univariate distribution.
"""
struct UnivariateFill{S<:Sampleable{Univariate}}
    sampleable::S
    dims::Dims{2}
end
function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:UnivariateFill})
    return rand(rng, d[].sampleable, d[].dims)
end
function Random.rand!(rng::AbstractRNG, a::AbstractMatrix, d::SamplerTrivial{<:UnivariateFill})
    return Random.rand!(rng, d[].sampleable, a)
end
Base.eltype(d::UnivariateFill) = Matrix{eltype(d.sampleable)}
Base.size(d::UnivariateFill) = d.dims
Base.getindex(d::UnivariateFill, _) = rand(d)

"""
    UnivariateArray{S<:Sampleable{Univariate}}

Fill each degree of freedom from a different Univariate distribution.

The size of the matrix of sampleables should match the system size.
"""
struct UnivariateArray{S<:Sampleable{Univariate}}
    sampleable::Matrix{S}
end
function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:UnivariateArray})
    out = similar(d[].sampleable, eltype(eltype(d[].sampleable)))
    return Random.rand!(rng, out, d)
end
function Random.rand!(rng::AbstractRNG, a::AbstractMatrix, d::SamplerTrivial{<:UnivariateArray})
    for I in eachindex(a, d[].sampleable)
        a[I] = rand(rng, d[].sampleable[I])
    end
    return a
end
Base.eltype(d::UnivariateArray) = Matrix{eltype(eltype(d.sampleable))}
Base.size(d::UnivariateArray) = size(d.sampleable)
Base.getindex(d::UnivariateArray, _) = rand(d)

"""
    FixedArray{S<:AbstractArray}

Return the same configuration every time.
"""
struct FixedArray{S<:AbstractArray}
    value::S
end
function Random.rand(::AbstractRNG, d::SamplerTrivial{<:FixedArray})
    return copy(d[].value)
end
function Random.rand!(::AbstractRNG, a::AbstractArray, d::SamplerTrivial{<:FixedArray})
    copy!(a, d[].value)
    return a
end
Base.eltype(d::FixedArray) = typeof(d.value)
Base.size(d::FixedArray) = size(d.value)
Base.getindex(d::FixedArray, _) = rand(d)

"""
    FixedValue{S<:Real}

Fill all degrees of freedom with the same value every time.
"""
struct FixedValue{S<:Real}
    value::S
    dims::Dims{2}
end
function Random.rand(::AbstractRNG, d::SamplerTrivial{<:FixedValue})
    return fill(d[].value, d[].dims)
end
function Random.rand!(::AbstractRNG, a::AbstractMatrix, d::SamplerTrivial{<:FixedValue})
    fill!(a, d[].value)
    return a
end
Base.eltype(d::FixedValue) = Matrix{typeof(d.value)}
Base.size(d::FixedValue) = d.dims
Base.getindex(d::FixedValue, _) = rand(d)

"""
    ConfigurationVector{S<:AbstractVector}

Sample from a provided vector of configurations.
"""
struct ConfigurationVector{S<:AbstractVector}
    configurations::S
end
function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:ConfigurationVector})
    return copy(rand(rng, d[].configurations))
end
function Random.rand!(rng::AbstractRNG, a::AbstractArray, d::SamplerTrivial{<:ConfigurationVector})
    copy!(a, rand(rng, d[].configurations))
    return a
end
Base.eltype(d::ConfigurationVector) = eltype(d.configurations)
Base.size(d::ConfigurationVector) = size(first(d.configurations))
Base.getindex(d::ConfigurationVector, i) = copy(d.configurations[i])

"""
    RingPolymerWrapper{S}

Wrap other distributions to convert them to ring polymer distributions.
"""
struct RingPolymerWrapper{S}
    sampleable::S
    dims::Dims{3}
    classical::Vector{Int}
end
function RingPolymerWrapper(sampleable, nbeads; classical=Int[])
    return RingPolymerWrapper(sampleable, (size(sampleable)..., nbeads), classical)
end
function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:RingPolymerWrapper})
    out = RingPolymerArray{eltype(eltype(d[].sampleable))}(undef, d[].dims; d[].classical)
    Random.rand!(rng, out, d)
    return out
end
function Random.rand!(rng::AbstractRNG, a::RingPolymerArray, d::SamplerTrivial{<:RingPolymerWrapper})
    for bead in eachbead(a)
        Random.rand!(rng, bead, d[].sampleable)
    end
    return a
end
Base.eltype(d::RingPolymerWrapper) = RingPolymerArray{eltype(eltype(d.sampleable))}
Base.size(d::RingPolymerWrapper) = d.dims
function Base.getindex(d::RingPolymerWrapper, i)
    out = RingPolymerArray{eltype(eltype(d.sampleable))}(undef, d.dims; d.classical)
    for bead in eachbead(out)
        copy!(bead, d.sampleable[i])
    end
    return out
end

end
