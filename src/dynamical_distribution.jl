
struct DynamicalDistribution{V,R}
    velocity::V
    position::R
    function DynamicalDistribution(velocity, position)
        size(velocity) == size(position) || throw(
            DimensionMismatch(
                "`velocity` and `position` sample size does not match: \
                $(size(velocity)) != $(size(position))"
            )
        )
        return new{typeof(velocity),typeof(position)}(velocity, position)
    end
end

function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:DynamicalDistribution})
    i = rand(rng, 1:lastindex(d[]))
    return d[][i]
end

function Base.getindex(d::DynamicalDistribution, i)
    return ComponentVector(v=d.velocity[i], r=d.position[i])
end

function Base.lastindex(d::DynamicalDistribution)
    return max(lastindex(d.velocity), lastindex(d.position))
end

function Base.iterate(d::DynamicalDistribution, state=1)
    return state > lastindex(d) ? nothing : (d[state], state+1)
end
