module NQCDSizeDistributionExt

using NQCDynamics

"""
    VelocityBoltzmann(temperature, sim::NQCDynamics.AbstractSimulation; center = zeros(size(sim)))

Simulation-aware version of VelocityBoltzmann. 
This will generate a VelocityBoltzmann distribution with the correct size for the simulation using the masses specified in the simulation.
"""
function VelocityBoltzmann(temperature, sim::NQCDynamics.AbstractSimulation; center = zeros(size(sim)))
	return VelocityBoltzmann(temperature, sim.atoms.masses, size(sim); center=center)
end

end # module