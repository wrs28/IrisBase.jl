"""
    isCartesian(sim)
"""
function CoordinateSystems.isCartesian(sim::Simulation)
    return isCartesian(sim.dis)
end
function CoordinateSystems.isCartesian(dis::Discretization)
    return isCartesian(dis.coordinate_system)
end


"""
    isPolar(sim)
"""
function CoordinateSystems.isPolar(sim::Simulation)
    return isPolar(sim.dis)
end
function CoordinateSystems.isPolar(dis::Discretization)
    return isPolar(dis.coordinate_system)
end
