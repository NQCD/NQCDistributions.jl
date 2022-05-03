
function VelocityBoltzmann(temperature, masses::AbstractVector, dims::Dims{2})
    return UnivariateArray([VelocityBoltzmann(temperature, masses[I[2]]) for I in CartesianIndices(dims)])
end
VelocityBoltzmann(temperature, mass) = Normal(0, sqrt(austrip(temperature)/mass))
