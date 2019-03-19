module CoordinateSystems

export CoordinateSystem,
Polar,
Cartesian,
isPolar,
isCartesian

abstract type CoordinateSystem end
struct Polar<:CoordinateSystem end
struct Cartesian<:CoordinateSystem end

"""
    bool = isCartesian(coordinate_system)
"""
isCartesian(coordinate_system::Cartesian) = true
isCartesian(coordinate_system) = false


"""
    bool = isPolar(coordinate_system)
"""
isPolar(coordinate_system::Polar) = true
isPolar(coordinate_system) = false

end # module
