
"""
    VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2})

Generate a Boltzmann of velocities for each degree of freedom.

# Arguments

* `temperature` - Atomic units or Unitful
* `masses` - Vector of masses for each atom
* `dims` - (ndofs, natoms). `natoms` must equal `length(masses)`
"""
function VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2})
    return UnivariateArray([VelocityBoltzmann(temperature, masses[I[2]]) for I in CartesianIndices(dims)])
end
VelocityBoltzmann(temperature, mass) = Normal(0, sqrt(austrip(temperature)/mass))
