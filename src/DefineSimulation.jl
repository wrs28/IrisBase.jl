module DefineSimulation

const DEFAULT_SUBSAMPLE_NUMBER = 5 # default number of samples used in sub-pixel smoothing (total number is square of this)

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
Statistics,
StaticArrays

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

include("DefineSimulation/domain_types.jl")
include("DefineSimulation/main_structs.jl")
include("DefineSimulation/construction.jl")

include("DefineSimulation/overloading.jl")
include("DefineSimulation/coordinates.jl")
include("DefineSimulation/bravais.jl")

include("DefineSimulation/plot_defaults.jl")
include("DefineSimulation/plotting.jl")

end # module
