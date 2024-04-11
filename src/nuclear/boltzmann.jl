
"""
    VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2};  center::AbstractVector = zeros(dims))

Generate a Boltzmann distribution of velocities for each degree of freedom.

# Arguments

* `temperature` - Atomic units or Unitful
* `masses` - Vector of masses for each atom
* `dims` - (ndofs, natoms). `natoms` must equal `length(masses)`
* `center` - Vector of dimension ndofs that provides centre of mass velocity offset
"""

function VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2};  center::Array{Float64} = zeros(dims))
    return UnivariateArray([VelocityBoltzmann(temperature, masses[I[2]]; center = center[I]) for I in CartesianIndices(dims)])
end
    
"""
    VelocityBoltzmann(temperature, mass; center = 0)
"""
function VelocityBoltzmann(temperature::Number, mass::Float64; center::Float64 = 0)
    return Normal(center, sqrt(austrip(temperature)/mass))
end
