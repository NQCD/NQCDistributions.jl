
"""
    VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2})

Generate a Boltzmann distribution of velocities for each degree of freedom.

# Arguments

* `temperature` - Atomic units or Unitful
* `masses` - Vector of masses for each atom
* `dims` - (ndofs, natoms). `natoms` must equal `length(masses)`
* `centre` - Vector of dimension ndofs that provides centre of mass velocity offset
"""
function VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2};  centre = nothing)
    if centre == nothing
        return UnivariateArray([VelocityBoltzmann(temperature, masses[I[2]]) for I in CartesianIndices(dims)])
    else
        if size(centre,1) ≠ dims[1]
            throw(error("centre needs to have length `ndof`"))
        end
        return UnivariateArray([VelocityBoltzmann(temperature, masses[I[2]]; centre = centre[I[1]]) for I in CartesianIndices(dims)])
    end
end
VelocityBoltzmann(temperature, mass; centre = 0) = Normal(centre, sqrt(austrip(temperature)/mass))
