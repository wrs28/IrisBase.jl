module DefineSimulation

using ..Bravais,
..CoordinateSystems,
..BoundaryConditions,
..DielectricFunctions,
..Dispersions,
..Shapes,
Formatting,
Interpolations,
IterTools,
RecipesBase,
Statistics

import ..BoundaryConditions: reorder, get_dim, get_side, apply_args
import ..DifferentialOperators: _oc_bls

export Domain,
System,
Discretization,
Boundary,
Channels,
Scattering,
TwoLevelSystem,
Simulation,
which_domain

# include("defaults.jl")
include("DefineSimulation/main_structs.jl")
include("DefineSimulation/construction.jl")

include("DefineSimulation/overloading.jl")
include("DefineSimulation/coordinates.jl")
include("DefineSimulation/bravais.jl")

include("DefineSimulation/iss.jl")

include("DefineSimulation/pretty_printing.jl")
include("DefineSimulation/plot_defaults.jl")
include("DefineSimulation/plotting.jl")

end # module
