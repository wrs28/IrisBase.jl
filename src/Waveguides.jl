# TODO: Test arbitrary profile waveguide, for now have only tested simple waveguide
module Waveguides

export add_halfspace_waveguide,
add_halfspace_waveguides,
add_planar_waveguide,
add_planar_waveguides,
add_pc_waveguide,
add_pc_waveguides

using ...Bravais,
...BoundaryConditions,
...Shapes,
...DefineSimulation,
..ConstructionToolsBase,
..PhotonicCrystal

include("Waveguides/waveguides_planar.jl")
include("Waveguides/waveguides_halfspace.jl")
include("Waveguides/waveguides_pc.jl")

include("DefineSimulation.jl/iss.jl")

end # module
