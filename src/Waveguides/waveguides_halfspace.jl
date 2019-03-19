"""
    waveguide_domains = planar_waveguides(;index, direction, waveguide_number)
"""
function halfspace_waveguide_domain(index::Number, direction::Symbol, waveguide_number::Int)
    if direction ∈ [:w, :l, :left, :west]
        which_asymptote = :left
        lattice = BravaisLattice(a=0)
    elseif direction ∈ [:e, :r, :right, :east]
        which_asymptote = :right
        lattice = BravaisLattice(a=0)
    elseif direction ∈ [:n, :u, :t, :north, :up, :top]
        which_asymptote = :top
        lattice = BravaisLattice(b=0)
    elseif direction ∈ [:s, :b, :d, :south, :bottom, :down]
        which_asymptote = :bottom
        lattice = BravaisLattice(b=0)
    else
        throw(ArgumentError("invalid waveguide direction $(direction)"))
    end

    wvg = build_domain(; which_asymptote=which_asymptote, which_waveguide=waveguide_number, domain_type=:halfspace_waveguide, domain_params=Dict(:n₁=>index), lattice=lattice)

    return wvg
end


"""
    sys = add_halfspace_waveguide(sys; index=1, direction)

    sim = add_halfspace_waveguide(sim; index=1, direction)
"""
function add_halfspace_waveguide(sys::System; index::Number=1, direction::Symbol, waveguide_number::Int=length(sys.waveguides)+1)
    waveguide_domains = halfspace_waveguide_domain(index, direction, waveguide_number)
    return System(vcat(waveguide_domains,sys.domains))
end
function add_halfspace_waveguide(sim::Simulation; index::Number=1, direction::Symbol)
    @assert 1 ∈ sim.dis.N "halfspace waveguides only for one-dimensional-systems"
    waveguide_number = length(sim.sys.waveguides)+1
    sys = add_halfspace_waveguide(sim.sys; index=index, direction=direction, waveguide_number=waveguide_number)
    return Simulation(sim; :sys => sys)
end


"""
    sim = add_halfspace_waveguide(sim, side::Symbol; index=1)
"""
function add_halfspace_waveguide(sim::Simulation, side::Symbol; index::Number=1)
    return add_halfspace_waveguide(sim; index=index, direction=side)
end


"""
    sim = add_halfspace_waveguides(sim; index=1)
"""
function add_halfspace_waveguides(sim::Simulation; index::Number=1)
    @assert 1 ∈ sim.dis.N "halfspace waveguides only for one-dimensional-systems"
    if sim.dis.N[1]==1
        sim = add_halfspace_waveguide(sim, :bottom; index=index)
        sim = add_halfspace_waveguide(sim, :top; index=index)
    else sim.dis.N[2]==1
        sim = add_halfspace_waveguide_left(sim, :left; index=index)
        sim = add_halfspace_waveguide_right(sim, :right; index=index)
    end
    return sim
end
