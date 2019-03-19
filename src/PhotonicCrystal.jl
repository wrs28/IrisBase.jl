module PhotonicCrystal

using ...Bravais,
...Shapes,
...DielectricFunctions,
...DefineSimulation,
...BoundaryConditions,
..ConstructionToolsBase

export build_pc_domain,
line_defect


"""
    sim = Simulation(domain; bnd=Boundary(bc=:p), dis, sct, tls, Na=1, Nb=1)

Simulation from a single domain.

Intended for case where `domain` has a finite lattice structure, in which case
    it defines the simulation to be `Na`×`Nb` unit cells.

Care must be taken with boundary layer specifications. It is recommended that
`bnd` be left unspecified.
"""
function DefineSimulation.Simulation(dom::Domain; bnd::Boundary=Boundary(bc=PeriodicBC), dis::Discretization, sct::Scattering=Scattering(), tls::TwoLevelSystem=TwoLevelSystem(), Na::Int=1, Nb::Int=1)
    lattice = dom.lattice

    # bnd = deepcopy(bnd)
    if !isinf(lattice.a)
        if typeof(bnd.bc[1][1])<:PeriodicBC
            ∂Ω1 = (lattice.x0,
                    lattice.x0+Na*lattice.a*(cos(lattice.α)-sin(lattice.α)*cot(lattice.β)))
        else
            ∂Ω1 = bnd.∂Ω[1]
        end
    end
    if !isinf(lattice.b)
        if typeof(bnd.bc[2][1])<:PeriodicBC
            ∂Ω2 = (lattice.y0,
                    lattice.y0+Nb*lattice.b*sin(lattice.β))
        else
            ∂Ω2 = bnd.∂Ω[2]
        end
    end
    return deepcopy(Simulation(System(dom), Boundary(bnd; ∂Ω=(∂Ω1,∂Ω2)), dis, tls, BravaisLattice(lattice; :a=>Na*lattice.a, :b=>Nb*lattice.b)))
end
"""
    sim = Simulation(domain, sim; Na=1, Nb=1)

Build simulation out of a single `domain`, with fields other than `bnd` the same
as in `sim`, except it is `Na`×`Nb` unit cells.

Intended for case where `domain` has a finite lattice structure.
"""
DefineSimulation.Simulation(dom::Domain, sim::Simulation; Na::Int=1, Nb::Int=1) = Simulation(dom; dis=sim.dis, tls=sim.tls, Na=Na, Nb=Nb)
"""
    sim = Simulation(domain_num, sim; Na=1, Nb=1)

Build simulation out of domain number `domain_num`, with fields other than `bnd`
the same as in `sim`, except it is `Na`×`Nb` unit cells.

Intended for case where specified domain has a finite lattice structure.
"""
Simulation(domain_num::Int, sim::Simulation; Na::Int=1, Nb::Int=1) = deepcopy(Simulation(sim.sys.domains[domain_num]; dis=sim.dis, tls=sim.tls, Na=Na, Nb=Nb))
"""
    domain = build_pc_domain
"""
function build_pc_domain(
            lattice::BravaisLattice,
            subdomains::Tsd = Subdomains();
            is_in_domain::Tsh = Universe(),
            domain_params::Td = Dict(:n₁ => 1, :n₂ => 0, :F => 0),
            domain_type::Symbol = :pc,
            domain_ε::Function = piecewise_constant_ε,
            domain_F::Function = piecewise_constant_F,
            which_asymptote::Symbol = :none,
            which_waveguide::Int = 0
            ) where Tsd<:Subdomains where Tsh<:AbstractShape where Td<:Dict

    return build_domain(
                subdomains;
                is_in_domain=is_in_domain,
                domain_params=domain_params,
                domain_type=domain_type,
                domain_ε=domain_ε,
                domain_F=domain_F,
                lattice=lattice,
                which_asymptote=which_asymptote,
                which_waveguide=which_waveguide
                )
end


##########################################################################################
### DEFECT CONSTRUCTION
##########################################################################################
"""
    line_defect(pc_domain, direction, N_start, N_stop, N_transverse, params, ε, F, subdomains; N_width=1)
"""
function line_defect(
    pc_domain::Tdom,
    direction::Symbol,
    N_start::Int,
    N_stop::Int,
    N_transverse::Int,
    params=pc_domain.domain_params,
    ε::Function=pc_domain.domain_ε,
    F::Function=pc_domain.domain_F,
    subdomains::Tsd=Subdomains();
    N_width::Int=1
    ) where Tdom<:Domain where Tsd<:Subdomains

    @assert direction ∈ [:a, :b] "invalid direction $(direction)."
    if direction == :a
        Na_start = N_start
        Na_stop = N_stop
        Nb_start = N_transverse
        Nb_stop = N_transverse+N_width-1
    else direction == :b
        Nb_start = N_start
        Nb_stop = N_stop
        Na_start = N_transverse
        Na_stop = N_transverse+N_width-1
    end
    lattice = pc_domain.lattice

    if Na_start==Na_stop && Nb_start==Nb_stop
        defect_type = :site_defect
    else
        defect_type = :line_defect
    end

    a, b = lattice.a*(Na_stop-Na_start+1), lattice.a*(Nb_stop-Nb_start+1)
    α, θ = lattice.β-lattice.α, lattice.α
    x0 = lattice.x0 + Na_start*lattice.a*lattice.v1[1] + lattice.x0 + Nb_start*lattice.b*lattice.v2[1]
    y0 = lattice.y0 + Na_start*lattice.a*lattice.v1[2] + lattice.x0 + Nb_start*lattice.b*lattice.v2[2]

    return build_domain(subdomains;
                is_in_domain = Parallelogram(a,b,α,x0,y0,θ),
                domain_params = params,
                domain_type = defect_type,
                domain_ε = ε,
                domain_F = F,
                lattice = lattice)
end


"""
    site_defect
"""
function site_defect(pc_domain::Domain, Na::Int, Nb::Int,
        params=pc_domain.domain_params, ε=pc_domain.domain_ε, F=pc_domain.domain_F,
        sudomains::Subdomains=Subdomains() )

    defect_domain = line_defect_domain(pc_domain, :a, Na, Na, Nb, params, ε, F, sudomains)
    return defect_domain
end
function site_defect(pc_domain::Domain, Na::Int, Nb::Array{Int},
    params=pc_domain.domain_params, ε=pc_domain.domain_ε, F=pc_domain.domain_F,
    sudomains::Subdomains=Subdomains() )

    return defect_domain.(Ref(pc_domain), Na, Nb, Ref(params), Ref(ε), Ref(F), Ref(sudomains))
end
function site_defect(pc_domain::Domain, Na::Array{Int}, Nb,
    params=pc_domain.domain_params, ε=pc_domain.domain_ε, F=pc_domain.domain_F,
    sudomains::Subdomains=Subdomains())

    return defect_domain.(Ref(pc_domain), Na, Nb, Ref(params), Ref(ε), Ref(F), Ref(sudomains))
end

#################################################################################################################
### SHAPE OVERLOADING
#################################################################################################################
"""
    Circle(R, lattice, referece, params[, subdomains; x0, y0, type, ε, F])
"""
function Shapes.Circle(R::Real, reference::Symbol,
            params::Dict,
            domain::Domain;
            x0::Real=0,
            y0::Real=0,
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    x0_new, y0_new = get_reference(reference, x0, y0, domain.lattice)
    for i ∈ eachindex(x0_new)
        domain = Circle(R, x0_new[i], y0_new[i], params, domain; type=type, ε=ε, F=F)
    end
    return domain
end


"""
    Ellipse(a, b, θ, lattice, referece, params[, subdomains; x0, y0, type, ε, F])
"""
function Shapes.Ellipse(a::Real, b::Real, θ::Real, reference::Symbol,
            params::Dict,
            domain::Domain;
            x0::Number=0,
            y0::Number=0,
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    x0_new, y0_new = get_reference(reference, x0, y0, domain.lattice)
    for i ∈ eachindex(x0_new)
        domain = Ellipse(a, b, x0_new[i], y0_new[i], θ, params, domain; type=type, ε=ε, F=F)
    end
    return domain
end


"""
    Square(a, θ, lattice, referece, params[, subdomains; x0, y0, type, ε, F])
"""
function Shapes.Square(a::Real, θ::Real, reference::Symbol,
            params::Dict,
            domain::Domain;
            x0::Number=0,
            y0::Number=0,
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    x0_new, y0_new = get_reference(reference, x0, y0, domain.lattice)
    for i ∈ eachindex(x0_new)
        domain = Square(a, x0_new[i], y0_new[i], θ, params, domain; type=type, ε=ε, F=F)
    end
    return domain
end


"""
    Rectangle(a, b, θ, lattice, referece, params[, subdomains; x0, y0, type, ε, F])
"""
function Shapes.Rectangle(a::Real, b::Real, θ::Real, reference::Symbol,
            params::Dict,
            domain::Domain;
            x0::Number=0,
            y0::Number=0,
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    x0_new, y0_new = get_reference(reference, x0, y0, domain.lattice)
    for i ∈ eachindex(x0_new)
        domain = Rectangle(a, b, x0_new[i], y0_new[i], θ, params, domain; type=type, ε=ε, F=F)
    end
    return domain
end


"""
    Parallelogram(a, b, α, θ, lattice, referece, params[, subdomains; x0, y0, type, ε, F])
"""
function Shapes.Parallelogram(a::Real, b::Real, α::Real, θ::Real, reference::Symbol,
            params::Dict,
            domain::Domain;
            x0::Number=0,
            y0::Number=0,
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    x0_new, y0_new = get_reference(reference, x0, y0, domain.lattice)
    for i ∈ eachindex(x0_new)
        domain = Parallelogram(a, b, α, x0_new[i], y0_new[i], θ, params, domain; type=type, ε=ε, F=F)
    end
    return domain
end


"""
    get_reference(reference::Symbol, x0, y0, lattice::BravaisLattice)
"""
function get_reference(reference::Symbol, x0::Real, y0::Real, lattice::BravaisLattice)
    x0_new, y0_new = [], []
    v1, v2 = lattice.v1, lattice.v2
    a, b = lattice.a, lattice.b
    if reference == :center
        push!(x0_new, x0 - v1[1]*a/2 + v2[1]*b/2)
        push!(y0_new, y0 - v1[2]*a/2 + v2[2]*b/2)

        push!(x0_new, x0 + v1[1]*a/2 - v2[1]*b/2)
        push!(y0_new, y0 + v1[2]*a/2 - v2[2]*b/2)

        push!(x0_new, x0 + 3v1[1]*a/2 + v2[1]*b/2)
        push!(y0_new, y0 + 3v1[2]*a/2 + v2[2]*b/2)

        push!(x0_new, x0 + v1[1]*a/2 + 3v2[1]*b/2)
        push!(y0_new, y0 + v1[2]*a/2 + 3v2[2]*b/2)

        push!(x0_new, x0 + v1[1]*a/2 + v2[1]*b/2)
        push!(y0_new, y0 + v1[2]*a/2 + v2[2]*b/2)
    elseif reference ∈ [:corner, :corners]
        push!(x0_new, v1[1]*a)
        push!(y0_new, v1[2]*a)

        push!(x0_new, v2[1]*b)
        push!(y0_new, v2[2]*b)

        push!(x0_new, v1[1]*a + v2[1]*b)
        push!(y0_new, v1[2]*a + v2[2]*b)

        push!(x0_new, 0)
        push!(y0_new, 0)
    elseif reference ∈ [:bottom_edge, :top_edge, :bottom, :top]
        push!(x0_new, v1[1]*a/2 + v2[1]*b)
        push!(y0_new, v1[2]*a/2 + v2[2]*b)

        push!(x0_new, v1[1]*a/2)
        push!(y0_new, v1[2]*a/2)
    elseif reference ∈ [:left_edge, :right_edge, :left, :right]
        push!(x0_new, v1[1]*a + v2[1]*b/2)
        push!(y0_new, v1[2]*a + v2[2]*b/2)

        push!(x0_new, v2[1]*b/2)
        push!(y0_new, v2[2]*b/2)
    elseif reference ∈ [:southwest, :sw]
        push!(x0_new, x0)
        push!(y0_new, y0)
    elseif reference ∈ [:south, :s]
        push!(x0_new, x0 + v1[1]*a/2)
        push!(y0_new, y0 + v1[1]*a/2)
    elseif reference ∈ [:southeast, :se]
        push!(x0_new, x0 + v1[1]*a)
        push!(y0_new, y0 + v1[2]*a)
    elseif reference ∈ [:east, :e]
        push!(x0_new, x0 + v1[1]*a + v2[1]*b/2)
        push!(y0_new, y0 + v1[2]*a + v2[2]*b/2)
    elseif reference ∈ [:northeast, :ne]
        push!(x0_new, x0 + v1[1]*a + v2[1]*b)
        push!(y0_new, y0 + v1[2]*a + v2[2]*b)
    elseif reference ∈ [:north, :n]
        push!(x0_new, x0 + v1[1]*a/2 + v2[1]*b)
        push!(y0_new, y0 + v1[2]*a/2 + v2[2]*b)
    elseif reference ∈ [:northwest, :nw]
        push!(x0_new, x0 + v2[1]*b)
        push!(y0_new, y0 + v2[2]*b)
    elseif reference ∈ [:west, :w]
        push!(x0_new, x0 + v2[1]*b/2)
        push!(y0_new, y0 + v2[2]*b/2)
    else
        throw(ArgumentError("invalid unit cell position reference $(reference).
        Should be one of :center, :corner or a compass direction"))
    end
    return x0_new, y0_new
end


end #module
