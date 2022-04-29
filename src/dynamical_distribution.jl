
struct DynamicalDistribution{V,R}
    velocity::V
    position::R
    rng::Xoshiro
    function DynamicalDistribution(velocity, position)
        size(velocity) == size(position) || throw(
            DimensionMismatch(
                "`velocity` and `position` sample size does not match: \
                $(size(velocity)) != $(size(position))"
            )
        )
        return new{typeof(velocity),typeof(position)}(velocity, position, Xoshiro())
    end
end

function DynamicalDistribution(velocity, position, dims)
    v = SampleableComponent(velocity, dims)
    r = SampleableComponent(position, dims)
    return DynamicalDistribution(v, r)
end

function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:DynamicalDistribution})

    if isindexable(d[].velocity) || isindexable(d[].position)
        i = rand(rng, 1:lastindex(d[]))
        return d[][i]
    else
        i = rand(rng, UInt)
        Random.seed!(d[].rng, i)
        return ComponentVector(v=rand(d[].rng, d[].velocity), r=rand(d[].rng, d[].position))
    end

end

function Base.getindex(d::DynamicalDistribution, i)

    if isindexable(d.velocity) && isindexable(d.position)
        return ComponentVector(v=d.velocity[i], r=d.position[i])

    elseif isindexable(d.velocity)
        Random.seed!(d.rng, i)
        return ComponentVector(v=d.velocity[i], r=rand(d.rng, d.position))

    elseif isindexable(d.position)
        Random.seed!(d.rng, i)
        return ComponentVector(v=rand(d.rng, d.velocity), r=d.position[i])

    else
        Random.seed!(d.rng, i)
        return ComponentVector(v=rand(d.rng, d.velocity), r=rand(d.rng, d.position))
    end

end

function Base.lastindex(d::DynamicalDistribution)
    if isindexable(d.velocity) && isindexable(d.position)
        lastindex(d.velocity) == lastindex(d.position) || throw(
            DimensionMismatch("Length of position and velocity distributions do not match.")
        )
        return lastindex(d.velocity)

    elseif isindexable(d.velocity)
        return lastindex(d.velocity)

    elseif isindexable(d.position)
        return lastindex(d.position)

    else
        return 1

    end
end
Base.length(d::DynamicalDistribution) = lastindex(d)

function Base.iterate(d::DynamicalDistribution, state=1)
    return state > lastindex(d) ? nothing : (d[state], state+1)
end