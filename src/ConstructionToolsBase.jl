module ConstructionToolsBase

using ...DielectricFunctions,
...Shapes,
...Bravais,
...DefineSimulation

export Subdomains,
add_subdomain,
build_domain

"""
    subdomains = Subdomains()
    subdomains = Subdomains(is_in_subdomain::Array{AbstractShape}, params::Array{Dict}, type::Array{AbstractShape}, ε::Array{Function}, F::Array{Function})

`Subdomains()` initializes empty region

`Subdomains(functions, params, types, ε's, F's)` creates region with specified parameters
"""
struct Subdomains{TSD,TP}
    is_in_subdomain::Array{TSD,1}
    params::Array{TP,1}
    type::Array{Symbol,1}
    ε::Array{Function,1}
    F::Array{Function,1}
    function Subdomains()
        new{AbstractShape,Dict{Symbol,Float64}}(AbstractShape[], Dict{Symbol,Float64}[], Symbol[], Function[], Function[])
    end
    function Subdomains(is_in_subdomain::Array{Tsd,1}, params::Array{Tp,1}, type, ε, F) where Tsd<:AbstractShape where Tp<:Dict
        new{Tsd,Tp}(is_in_subdomain, params, type, ε, F)
    end
end
Base.show(io::IO, sbd::Subdomains) = begin
    print(io,
    typeof(sbd), " with ", length(sbd.is_in_subdomain), " subdomains:")
    temp0 = [["\n\t\tshape: ", sbd.is_in_subdomain[i]] for i ∈ eachindex(sbd.params)]
    temp1 = [["\n\t\tε: ", sbd.ε[i]] for i ∈ eachindex(sbd.ε)]
    temp2 = [["\n\t\tF: ", sbd.F[i]] for i ∈ eachindex(sbd.F)]
    temp3 = [["\n\t\tparams: ", sbd.params[i]] for i ∈ eachindex(sbd.params)]
    temp4 = [["\n\t\ttype: ", sbd.type[i]] for i ∈ eachindex(sbd.type)]
    for i ∈ eachindex(temp1)
        print(io, "\n\tsubdomain ", i)
        print(io, temp0[i]...)
        print(io, temp1[i]...)
        print(io, temp2[i]...)
        print(io, temp3[i]...)
        print(io, temp4[i]...)
    end
end


################################################################################
# GENERIC REGION CONSTRUCTION
################################################################################
"""
    add_subdomain(is_in_subdomain::AbstractShape, params::Dict, type::Symbol, ε::Function, F::Function, old_subdomains::Subdomains)
"""
function add_subdomain(
                    is_in_subdomain::AbstractShape,
                    params = Dict(),
                    type::Symbol = :background,
                    ε::Function = piecewise_constant_ε,
                    F::Function = piecewise_constant_F,
                    old_subdomains::Subdomains = Subdomains()
                    )

    return Subdomains(
                vcat(is_in_subdomain, old_subdomains.is_in_subdomain),
                vcat(params, old_subdomains.params),
                vcat(type, old_subdomains.type),
                vcat(ε, old_subdomains.ε),
                vcat(F, old_subdomains.F)
                )
end
"""
    add_subdomain(domain::Domain; is_in_subdomain::AbstractShape, params::Dict, type::Symbol, ε::Function, F::Function)
"""
function add_subdomain(
            domain::Domain;
            is_in_subdomain::AbstractShape,
            params = Dict(),
            type::Symbol = :background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F
            )

    subdomains =  Subdomains(
                vcat(is_in_subdomain, domain.is_in_subdomain),
                vcat(params, domain.subdomain_params),
                vcat(type, domain.subdomain_type),
                vcat(ε, domain.subdomain_ε),
                vcat(F, domain.subdomain_F)
                )

    return build_domain(subdomains;
        is_in_domain = domain.is_in_domain,
        domain_params = domain.domain_params,
        domain_type = domain.domain_type,
        domain_ε = domain.domain_ε,
        domain_F = domain.domain_F,
        lattice = domain.lattice,
        which_asymptote = domain.which_asymptote,
        which_waveguide = domain.which_waveguide
        )
end
"""
    domain = add_subdomain(domain, subdomains)
    domain = add_subdomain(subdomains, domain)
"""
function add_subdomain(domain::Domain, subdomains::Subdomains)

    subdomains =  Subdomains(
                vcat(subdomains.is_in_subdomain, domain.is_in_subdomain),
                vcat(subdomains.params, domain.subdomain_params),
                vcat(subdomains.type, domain.subdomain_type),
                vcat(subdomains.ε, domain.subdomain_ε),
                vcat(subdomains.F, domain.subdomain_F)
                )

    return build_domain(subdomains;
        is_in_domain = domain.is_in_domain,
        domain_params = domain.domain_params,
        domain_type = domain.domain_type,
        domain_ε = domain.domain_ε,
        domain_F = domain.domain_F,
        lattice = domain.lattice,
        which_asymptote = domain.which_asymptote,
        which_waveguide = domain.which_waveguide
        )
end
function add_subdomain(subdomains::Subdomains, domain::Domain)
    return add_subdomaiun(domain, subdomain)
end



"""
    build_domain
"""
function build_domain(
            subomains::Tsd = Subdomains();
            is_in_domain::Tsh = Universe(),
            domain_params::Td = Dict(:n₁ => 1, :n₂ => 0, :F => 0),
            domain_type::Symbol = :background,
            domain_ε::Function = piecewise_constant_ε,
            domain_F::Function = piecewise_constant_F,
            lattice::BravaisLattice = BravaisLattice(),
            which_asymptote::Symbol = :none,
            which_waveguide::Int = 0
            ) where Tsd<:Subdomains where Tsh<:AbstractShape where Td<:Dict

    return Domain(; is_in_domain=is_in_domain,
                    domain_params=domain_params,
                    domain_type=domain_type,
                    domain_ε=domain_ε,
                    domain_F=domain_F,
                    is_in_subdomain = subomains.is_in_subdomain,
                    subdomain_params = subomains.params,
                    subdomain_type = subomains.type,
                    subdomain_ε = subomains.ε,
                    subdomain_F = subomains.F,
                    lattice=lattice,
                    which_asymptote=which_asymptote,
                    which_waveguide=which_waveguide
                    )
end

################################################################################
# SHAPE OVERLOADING
################################################################################
"""
    subdomains = Circle(R ,x0 ,y0 ,params::Dict, subdomains::Subdomains; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Circle(R::Number, x0::Number, y0::Number,
            params::Dict,
            subdomains::Subdomains;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return add_subdomain(Circle(R,x0,y0), params, type, ε, F, subdomains)
end
"""
    domain = Circle(R, x0, y0, params::Dict, domain::Domain; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Circle(R::Number, x0::Number, y0::Number,
            params::Dict,
            domain::Domain;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    subdomains = Circle(R, x0, y0, params, Subdomains(); type=type, ε=ε, F=F)
    return add_subdomain(domain, subdomains)
end
"""
    domain = Circle(R, x0, y0, params::Dict; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Circle(R::Number, x0::Number, y0::Number,
            params::Dict;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return build_domain(;is_in_domain=Circle(R,x0,y0), domain_params=params,
                            domain_type=type, domain_ε=ε, domain_F=F)
end


"""
    subdomains = Ellipse(a, b, x0, y0, θ, params::Dict, subdomains::Subdomains; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Ellipse(a::Number, b::Number, x0::Number, y0::Number, θ::Number,
            params::Dict,
            subdomains::Subdomains;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return add_subdomain(Ellipse(a,b,x0,y0,θ), params, type, ε, F, subdomains)
end
"""
    domain = Ellipse(a, b, x0, y0, θ, params::Dict, domain::Domain; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Ellipse(a::Number, b::Number, x0::Number, y0::Number, θ::Number,
            params::Dict,
            domain::Domain;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    subdomains = Ellipse(a, b, x0, y0, θ, params, Subdomains(); type=type, ε=ε, F=F)
    return add_subdomain(domain, subdomains)
end
"""
    domain = Ellipse(a, b, x0, y0, θ, params::Dict; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Ellipse(a::Number, b::Number, x0::Number, y0::Number, θ::Number,
            params::Dict;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return build_domain(;is_in_domain=Ellipse(a,b,x0,y0,θ), domain_params=params,
                            domain_type=type, domain_ε=ε, domain_F=F)
end


"""
    subdomains = Square(a, x0, y0, θ, params::Dict, subdomains::Subdomains; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Square(a::Number, x0::Number, y0::Number, θ::Number,
            params::Dict,
            subdomains::Subdomains;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return add_subdomain(Square(a,x0,y0,θ), params, type, ε, F, subdomains)
end
"""
    domain = Square(a, x0, y0, θ, params::Dict, domain::Domain; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Square(a::Number, x0::Number, y0::Number, θ::Number,
            params::Dict,
            domain::Domain;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    subdomains = Square(a, x0, y0, θ, params, Subdomains(); type=type, ε=ε, F=F)
    return add_subdomain(domain, subdomains)
end
"""
    domain = Square(a, x0, y0, θ, params::Dict; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Square(a::Number, x0::Number, y0::Number, θ::Number,
            params::Dict;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return build_domain(;is_in_domain=Square(a,x0,y0,θ), domain_params=params,
                            domain_type=type, domain_ε=ε, domain_F=F)
end


"""
    subdomains = Rectangle(a, b, x0, y0, θ, params::Dict, subdomains::Subdomains; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Rectangle(a::Number, b::Number, x0::Number, y0::Number, θ::Number,
            params::Dict,
            subdomains::Subdomains;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return add_subdomain(Rectangle(a, b, x0, y0, θ), params, type, ε, F, subdomains)
end
"""
    domain = Rectangle(a, b, x0, y0, θ, params::Dict, domain::Domain; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Rectangle(a::Number, b::Number, x0::Number, y0::Number, θ::Number,
            params::Dict,
            domain::Domain;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    subdomains = Rectangle(a, b, x0, y0, θ, params, Subdomains(); type=type, ε=ε, F=F)
    return add_subdomain(domain, subdomains)
end
"""
    domain = Rectangle(a, b, x0, y0, θ, params::Dict; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Rectangle(a::Number, b::Number, x0::Number, y0::Number, θ::Number,
            params::Dict;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return build_domain(;is_in_domain=Rectangle(a,b,x0,y0,θ), domain_params=params,
                            domain_type=type, domain_ε=ε, domain_F=F)
end


"""
    subdomains = Parallelogram(a, b, α, x0, y0, θ, params::Dict, subdomains::Subdomains; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Parallelogram(a::Number, b::Number, α::Number, x0::Number, y0::Number, θ::Number,
            params::Dict,
            subdomains::Subdomains;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return add_subdomain(Parallelogram(a, b, α, x0, y0, θ), params, type, ε, F, subdomains)
end
"""
    domain = Parallelogram(a, b, α, x0, y0, θ, params::Dict, domain::Domain; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Parallelogram(a::Number, b::Number, α::Number, x0::Number, y0::Number, θ::Number,
            params::Dict,
            domain::Domain;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    subdomains = Parallelogram(a, b, α, x0, y0, θ, params, Subdomains(); type=type, ε=ε, F=F)
    return add_subdomain(domain, subdomains)
end
"""
    domain = Parallelogram(a, b, α, x0, y0, θ, params::Dict; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Parallelogram(a::Number, b::Number, α::Number, x0::Number, y0::Number, θ::Number,
            params::Dict;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return build_domain(;is_in_domain=Parallelogram(a,b,α,x0,y0,θ), domain_params=params,
                            domain_type=type, domain_ε=ε, domain_F=F)
end


"""
    subdomains = Universe(params::Dict, subdomains::Subdomains; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Universe(
            params::Dict,
            subdomains::Subdomains;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return add_subdomain(Universe(), params, type, ε, F, subdomains)
end
"""
    domain = Universe(params::Dict; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Universe(
            params::Dict;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return build_domain(;is_in_domain=Universe(), domain_params=params,
                            domain_type=type, domain_ε=ε, domain_F=F)
end
"""
    domain = Universe(params::Dict; type=:background, ε=piecewise_constant_ε, F=piecewise_constant_F)
"""
function Shapes.Universe(
            n1::Number;
            type::Symbol=:background,
            ε::Function = piecewise_constant_ε,
            F::Function = piecewise_constant_F)

    return build_domain(;is_in_domain=Universe(), domain_params=Dict(:n1=>n1),
                            domain_type=type, domain_ε=ε, domain_F=F)
end

end #module
