"""
    waveguide_domains = planar_waveguides(;width, index, direction, x0, y0, waveguide_number)

returns vector of Domains
"""
function planar_waveguide_domains(width::Number, index::Number, direction::Symbol, x0::Number, y0::Number, waveguide_number::Int, background_index::Number)
    if direction ∈ [:ew, :lr, :h, :horizontal, :x]
        bulk_domain_type = :bulk_planar_waveguide_x
        which_asymptote1 = :left
        which_asymptote2 = :right
        a = 1e5
        domain = Rectangle(a=a, b=width, x0=x0-a/2, y0=y0-width/2, Dict(:n₁=>index))
    elseif direction ∈ [:ns, :ud, :tb, :v, :vertical, :y]
        bulk_domain_type = :bulk_planar_waveguide_y
        which_asymptote1 = :top
        which_asymptote2 = :bottom
        b=1e5
        domain = Rectangle(a=width, b=b, x0=x0-width/2, y0=y0-b/2, Dict(:n₁=>index))
    else
        throw(ArgumentError("invalid waveguide direction $(direction)"))
    end

    bgd1 = Domain(; :domain_type => :planar_waveguide_background, :n₁ => background_index, :which_waveguide => waveguide_number, :which_asymptote => which_asymptote1)
    bgd2 = Domain(; :domain_type => :planar_waveguide_background, :n₁ => background_index, :which_waveguide => waveguide_number+1, :which_asymptote => which_asymptote2)
    wvg1 = Domain(domain; :domain_type => :planar_waveguide, :which_waveguide => waveguide_number, :which_asymptote => which_asymptote1)
    wvg2 = Domain(domain; :domain_type => :planar_waveguide, :which_waveguide => waveguide_number+1, :which_asymptote => which_asymptote2)
    wvg_bulk = Domain(domain; :domain_type => bulk_domain_type)

    return vcat(bgd1, bgd2, wvg1, wvg2, wvg_bulk)
end


"""
    new_sys = add_planar_waveguide(sys; width, index, direction, position)

    new_sim = add_planar_waveguide(sim; width, index, direction, position)
"""
function add_planar_waveguide(sys::System; width::Number, index::Number, direction::Symbol, x0::Number, y0::Number, background_index::Number=1)
    waveguide_number = length(sys.waveguides)+1
    waveguide_domains = planar_waveguide_domains(width, index, direction, x0, y0, waveguide_number, background_index)
    return System(vcat(waveguide_domains,sys.domains))
end
function add_planar_waveguide(sim::Simulation; width::Number, index::Number, direction::Symbol, x0::Number, y0::Number, background_index::Number=sim.sys.domains[end].params[:n₁])
    sys = add_planar_waveguide(sim.sys; width=width, index=index, direction=direction, x0=x0, y0=y0, background_index=background_index)
    return Simulation(sim; :sys => sys)
end


"""
    new_sim = add_planar_waveguides(sim, dim; width, index, position)
"""
function add_planar_waveguides(sim::Simulation, dim::Int; width::Number, index::Number, direction::Symbol, x0::Number, y0::Number, background_index::Number=sim.sys.domains[end].params[:n₁])
    if dim == 1
        return add_horizontal_planar_waveguides(sim; width=width, index=index, y0=position, background_index=background_index)
    elseif dim == 2
        return add_vertical_planar_waveguides(sim; width=width, index=index, x0=position, background_index=background_index)
    else
        throw(ArgumentError("invalid dimension $dim"))
    end
end
function add_planar_waveguides(sim::Simulation, dir::Symbol; width::Number, index::Number, direction::Symbol, x0::Number, y0::Number, background_index::Number=sim.sys.domains[end].params[:n₁])
    @assert dir ∈ [:v, :V, :vertical, :Vertical, :vert, :Vert, :h, :H, :horizontal, :Horizontal, :horiz, :Horiz] "unrecognized direction $direction"
    if dir ∈ [:v, :V, :vertical, :Vertical, :vert, :Vert]
        return add_planar_waveguides(sim, 2; width=width, index=index, direction=direction, x0=x0, y0=y0, background_index=background_index)
    else
        return add_planar_waveguides(sim, 1; width=width, index=index, direction=direction, x0=x0, y0=y0, background_index=background_index)
    end
end


"""
    new_sim = add_horizontal_planar_waveguides(sim; width, index, y0)

    new_sys = add_horizontal_planar_waveguides(sys; width, index, y0)
"""
function add_horizontal_planar_waveguides(sys::System; width, index, y0, background_index=1)
    sys = add_planar_waveguide(sys; width=width, index=index, direction=:h, x0=0, y0=y0, background_index=background_index)
    return sys
end
function add_horizontal_planar_waveguides(sim::Simulation; width, index, y0, background_index=sim.sys.domains[end][:n₁])
    sys = add_horizontal_planar_waveguides(sim.sys, width=width, index=index, y0=y0, background_index=background_index)
    return Simulation(sim; :sys => sys)
end


"""
    new_sim = add_vertical_planar_waveguides(sim; width, index, x0)

    new_sys = add_vertical_planar_waveguides(sys; width, index, x0)
"""
function add_vertical_planar_waveguides(sys::System; width::Number, index::Number, x0::Number, background_index::Number=1)
    sys = add_planar_waveguide(sys; width=width, index=index, direction=:v, x0=x0, y0=0, background_index=background_index)
    return sys
end
function add_vertical_planar_waveguides(sim::Simulation; width::Number, index::Number, x0::Number, background_index::Number=sim.sys.domains[end][:n₁])
    sys = add_vertical_planar_waveguides(sim.sys; width=width, index=index, x0=x0, background_index=background_index)
    return Simulation(sim; :sys => sys)
end


"""
    new_sim = add_planar_waveguide(sim, side::Symbol; index=1)

    new_sys = add_planar_waveguide(sys, side::Symbol; index=1)
"""
function add_planar_waveguide(sim::Simulation, side::Symbol; index::Number=1)
    throw(ErrorException("this is intended for metallic waveguides, haven't done yet"))
    if side == :top
        direction = :n
        x0 = sim.bnd.∂Ω[1,1]
        y0 = sim.bnd.∂Ω[2,2]
        width = sim.bnd.∂Ω[2,1]-sim.bnd.∂Ω[1,1]
    elseif side == :bottom
        direction = :s
        x0 = sim.bnd.∂Ω[1,1]
        y0 = sim.bnd.∂Ω[1,2]
        width = sim.bnd.∂Ω[2,1]-sim.bnd.∂Ω[1,1]
    elseif side == :left
        direction = :w
        x0 = sim.bnd.∂Ω[1,1]
        y0 = sim.bnd.∂Ω[1,2]
        width = sim.bnd.∂Ω[2,2]-sim.bnd.∂Ω[1,2]
    elseif side == :right
        direction = :e
        x0 = sim.bnd.∂Ω[2,1]
        y0 = sim.bnd.∂Ω[1,2]
        width = sim.bnd.∂Ω[2,2]-sim.bnd.∂Ω[1,2]
    else
        throw(ArgumentError("invalide side $side, must be one of :top, :bottom, :left, :right"))
    end
    return add_planar_waveguide(simsys; width=width, index=index, direction=direction, x0=x0, y0=y0)
end
