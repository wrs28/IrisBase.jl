module ConstructionTools

export Subdomains,
add_subdomain,
build_domain,
build_pc_domain,
line_defect,
site_defect,
add_halfspace_waveguide,
add_halfspace_waveguides,
add_planar_waveguide,
add_planar_waveguides,
add_pc_waveguide,
add_pc_waveguides

include("ConstructionToolsBase.jl")
using .ConstructionToolsBase
#
# include("PhotonicCrystal.jl")
# using .PhotonicCrystal
#
# include("Waveguides.jl")
# using .Waveguides

end #module
