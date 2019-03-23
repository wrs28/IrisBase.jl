"""
    lattice = BravaisLattice(sim::Simulation)

Rectangular lattice with same size as defined in `sim.bnd` along periodic directions.
"""
# Bravais.BravaisLattice(sim::Simulation) = BravaisLattice(sim.bnd)


"""
    lattice = BravaisLattice(bnd::Boundary)

Rectangular lattice with same size as defined in `bnd` along periodic directions.
"""
function Bravais.BravaisLattice(bnd::Boundary)
    if typeof(bnd.bc[1][1])<:FloquetBC && typeof(bnd.bc[1][2])<:FloquetBC
        a = bnd.∂Ω[1][2]-bnd.∂Ω[1][1]
    elseif (typeof(bnd.bc[1][1])<:FloquetBC && !(typeof(bnd.bc[1][2])<:FloquetBC)) || (typeof(bnd.bc[1][2])<:FloquetBC && !(typeof(bnd.bc[1][1])<:FloquetBC))
        throw(ArgumentError("only one boundary condition along dimension 1 is periodic."))
    else
        a = Inf
    end
    if typeof(bnd.bc[2][1])<:FloquetBC && typeof(bnd.bc[2][2])<:FloquetBC
        b = bnd.∂Ω[2][2]-bnd.∂Ω[2][1]
    elseif (typeof(bnd.bc[2][1])<:FloquetBC && !(typeof(bnd.bc[2][2])<:FloquetBC)) || (typeof(bnd.bc[2][2])<:FloquetBC && !(typeof(bnd.bc[2][1])<:FloquetBC))
        throw(ArgumentError("only one boundary condition along dimension 2 is periodic."))
    else
        b = Inf
    end
    return BravaisLattice(a=a, b=b)
end


"""
    xb, yb = bravais_coordinates_unit_cell(x, y, domain_index, sys)

maps cartesian (`x`,`y`) into cartesian (`xb`,`yb`) unit cell specified by
sys.domains[domain_index].lattice
"""
# Bravais.bravais_coordinates_unit_cell(x, y, domain::Int, system::System) =
    # bravais_coordinates_unit_cell(x, y, system.domains[domain])


"""
    xb, yb = bravais_coordinates_unit_cell(x, y, domain)

maps cartesian (`x`,`y`) into cartesian (`xb`,`yb`) unit cell specified by
domain.lattice
"""
# Bravais.bravais_coordinates_unit_cell(x, y, domain::Domain) =
#     bravais_coordinates_unit_cell(x, y, domain.lattice)
# function Bravais.bravais_coordinates_unit_cell(x, y, domains::Array)
#     lattices = [domains[i].lattice for i ∈ eachindex(domains)]
#     #lattices = Array{BravaisLattice}(undef,size(x,1),size(y,2))
#     # for i ∈ CartesianIndices(domains)
#         # lattices[i] = domains[i].lattice
#     # end
#     return bravais_coordinates_unit_cell(x, y, lattices)
# end


"""
    xb, yb = bravais_coordinates_unit_cell(x, y, sim)

maps cartesian (`x`,`y`) into cartesian (`xb`,`yb`) unit cell specified by
sim.lat
"""
# Bravais.bravais_coordinates_unit_cell(x, y, sim::Simulation) =
#     bravais_coordinates_unit_cell(x,y,sim.lat)



"""
    bravais_coordinates_unit_cell!(xb, yb, x, y, domain_index, sys)

maps cartesian (`x`,`y`) into pre-allocated cartesian (`xb`,`yb`) unit cell specified by
sys.domains[domain_index].lattice
"""
# Bravais.bravais_coordinates_unit_cell!(xb::Array, yb::Array, x::AbstractArray, y::AbstractArray, domains::Array, system::System) =
    # bravais_coordinates_unit_cell!(xb, yb, x, y, system.domains[domains])


"""
    bravais_coordinates_unit_cell(xb, yb, x, y, domain)

maps cartesian (`x`,`y`) into pre-allocated cartesian (`xb`,`yb`) unit cell specified by
domain.lattice
"""
function Bravais.bravais_coordinates_unit_cell!(xb::AbstractArray, yb::AbstractArray, x::AbstractArray, y::AbstractArray, domains::AbstractArray, sys::System)
    for i ∈ eachindex(xb)
        xb[i],yb[i] = bravais_coordinates_unit_cell(x[i],y[i],sys.domains[domains[i]].lattice)
    end
    return nothing
end


"""
    bravais_coordinates_unit_cell!(xb, yb, x, y, sim)

maps cartesian (`x`,`y`) into pre-allocated cartesian (`xb`,`yb`) unit cell specified by
sim.lat
"""
# Bravais.bravais_coordinates_unit_cell!(xb, yb, x, y, sim::Simulation) =
#     bravais_coordinates_unit_cell!(xb, yb, x,y,sim.lat)


"""
    p1, p2 = bravais_coordinates(x, y, domain_index, sys)

coordinates in bravais frame specified by sys.domains[domain_index].lattice
(i.e. (x,y) = p1*v1 + p2*v2)
"""
# Bravais.bravais_coordinates(x, y, domain, system::System) =
#     bravais_coordinates(x, y, system.domains[domain])
"""
    p1, p2 = bravais_coordinates(x, y, domain)

coordinates in bravais frame specified by domain.lattice
(i.e. (x,y) = p1*v1 + p2*v2)
"""
# Bravais.bravais_coordinates(x, y, domain::Domain) =
#     bravais_coordinates(x, y, domain.lattice)

"""
    p1, p2 = bravais_coordinates(x, y, sim)

coordinates in bravais frame specified by sim.lat
(i.e. (x,y) = p1*v1 + p2*v2)
"""
# Bravais.bravais_coordinates(x, y, sim::Simulation) =
#     bravais_coordinates(x, y, sim.lat)
