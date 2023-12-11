
"""
    VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2};  centre = nothing)

Generate a Boltzmann distribution of velocities for each degree of freedom.

# Arguments

* `temperature` - Atomic units or Unitful
* `masses` - Vector of masses for each atom
* `dims` - (ndofs, natoms). `natoms` must equal `length(masses)`
* `centre` - Vector of dimension ndofs that provides centre of mass velocity offset
"""
function VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2};  centre = zeros(dims))
    return UnivariateArray([VelocityBoltzmann(temperature, masses[I[2]]; centre = centre[I]) for I in CartesianIndices(dims)])
end
VelocityBoltzmann(temperature, mass; centre = 0) = Normal(centre, sqrt(austrip(temperature)/mass))

"""
    +(dist::VelocityBoltzmann, translation::AbstractArray)

Add a permanent translation to a Velocity distribution, e.g. for a centre of mass translation component. 
"""
function Base.:+(dist::UnivariateArray, translation::AbstractArray)
    if size(dist)!=size(translation)
        throw(DimensionMismatch("Distribution and translation array must have the same size."))
    end
    return UnivariateArray([Normal(dist.sampleable[I].μ+translation[I], dist.sampleable[I].σ) for I in CartesianIndices(dist.sampleable)])
end

"""
    +(dist::Normal, translation::Number)

Add a permanent translation to a Velocity distribution. 
"""
function Base.:+(dist::Normal, translation::Number)
    return Normal(dist.μ+convert(Float64, translation), dist.σ)
end
