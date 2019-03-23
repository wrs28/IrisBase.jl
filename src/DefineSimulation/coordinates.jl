"""
    isCartesian(sim)
"""
CoordinateSystems.isCartesian(sim::Simulation) = isCartesian(sim.dis)
CoordinateSystems.isCartesian(dis::Discretization) = isCartesian(dis.coordinate_system)


"""
    isPolar(sim)
"""
CoordinateSystems.isPolar(sim::Simulation) = isPolar(sim.dis)
CoordinateSystems.isPolar(dis::Discretization) = isPolar(dis.coordinate_system)
