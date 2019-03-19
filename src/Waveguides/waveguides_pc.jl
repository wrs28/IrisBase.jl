"""
    waveguide_domains = defect_waveguide(domain::Domain, x0, y0, direction, width)

domain should be of the :pc type.
"""
function pc_waveguide_domains(pc_domain; width, direction, x0=0, y0=0, waveguide_number)

        xb0, yb0 = bravais_coordinates(x0, y0, pc_domain)
        if direction ∈ [:u, :up, :n, :north, :t, :top]
            axis = :b
            N1_bulk = floor(Int,yb0/pc_domain.lattice.b)
            N2_bulk = +Int(1e4)
            Nt = floor(Int,xb0/pc_domain.lattice.a)
            N_width = ceil(Int,width/pc_domain.lattice.a)
            bulk_domain_type = :bulk_pc_waveguide_y
            which_asymptote = :top
        elseif direction ∈ [:b, :bottom, :s, :south, :d, :down]
            axis = :b
            N1_bulk = -Int(1e4)
            N2_bulk = floor(Int,yb0/pc_domain.lattice.b)
            Nt = floor(Int,xb0/pc_domain.lattice.a)
            N_width = ceil(Int,width/pc_domain.lattice.a)
            bulk_domain_type = :bulk_pc_waveguide_y
            which_asymptote = :bottom
        elseif direction ∈ [:r, :right, :e, :east]
            axis = :a
            N1_bulk = floor(Int,xb0/pc_domain.lattice.a)
            N2_bulk = +Int(1e4)
            Nt = floor(Int,yb0/pc_domain.lattice.b)
            N_width = ceil(Int,width/pc_domain.lattice.b)
            bulk_domain_type = :bulk_pc_waveguide_x
            which_asymptote = :right
        elseif direction ∈ [:l, :left, :w, :west]
            axis = :a
            N1_bulk = -Int(1e4)
            N2_bulk = floor(Int,xb0/pc_domain.lattice.a)
            Nt = floor(Int,yb0/pc_domain.lattice.b)
            N_width = ceil(Int,width/pc_domain.lattice.b)
            bulk_domain_type = :bulk_pc_waveguide_x
            which_asymptote = :left
        elseif direction ∈ [:h, :horizontal, :x]
            axis = :a
            N1_bulk = -Int(1e4)
            N2_bulk = +Int(1e4)
            Nt = floor(Int,yb0/pc_domain.lattice.b)
            N_width = ceil(Int,width/pc_domain.lattice.b)
            bulk_domain_type = :bulk_pc_waveguide_x
            which_asymptote1 = :left
            which_asymptote2 = :right
        elseif direction ∈ [:v, :vertical, :y]
            axis = :b
            N1_bulk = -Int(1e4)
            N2_bulk = +Int(1e4)
            Nt = floor(Int,xb0/pc_domain.lattice.a)
            N_width = ceil(Int,width/pc_domain.lattice.a)
            bulk_domain_type = :bulk_pc_waveguide_y
            which_asymptote1 = :bottom
            which_asymptote2 = :top
        else
            throw(ArgumentError("invalid waveguide direction $(direction)"))
        end

        defect_bulk = line_defect(pc_domain, axis, N1_bulk, N2_bulk, Nt; N_width=N_width)
        wvg_bulk = vcat(Domain(defect_bulk; :domain_type => bulk_domain_type))

        defect_bnd = line_defect(pc_domain, axis, -Int(1e4), Int(1e4), Nt; N_width=N_width)

        if direction ∈ [:v, :vertical, :y, :x, :horizontal, :h]
            wvg = vcat(
                    Domain(defect_bnd; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote1, :domain_type => :pc_waveguide),
                    Domain(pc_domain; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote1, :domain_type => :pc_waveguide_background, :is_in_domain => Universe()),
                    Domain(defect_bnd; :which_waveguide => waveguide_number+1, :which_asymptote => which_asymptote2, :domain_type => :pc_waveguide),
                    Domain(pc_domain; :which_waveguide => waveguide_number+1, :which_asymptote => which_asymptote2, :domain_type => :pc_waveguide_background, :is_in_domain => Universe())
                    )
        else
            wvg = vcat(
                    Domain(defect_bnd; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote, :domain_type => :pc_waveguide),
                    Domain(pc_domain; :which_waveguide => waveguide_number, :which_asymptote => which_asymptote, :domain_type => :pc_waveguide_background, :is_in_domain => Universe())
                    )

        end

    return vcat(wvg, wvg_bulk)
end


"""
    sim = add_defect_waveguide(sim::Simulation; x0, y0, direction, width)
    sys = add_defect_waveguide(sys::System; x0, y0, direction, width)
"""
function add_pc_waveguide(sys::System; width::Number, direction::Symbol, x0::Number, y0::Number)
    waveguide_number = length(sys.waveguides)+1
    temp_sys = System(sys.domains[.!isWaveguide.([sys.domains[i] for i ∈ eachindex(sys.domains)]) .& .!isBulkWaveguide.([sys.domains[i] for i ∈ eachindex(sys.domains)]) .& .!isDefect.([sys.domains[i] for i ∈ eachindex(sys.domains)]) ])
    domain = which_domain(x0,y0, Boundary(∂Ω=[-Inf Inf;-Inf Inf],bc=DirichletBC), temp_sys)
    waveguide_domains = pc_waveguide_domains(temp_sys.domains[domain]; x0=x0, y0=y0, direction=direction, width=width, waveguide_number=waveguide_number)
    return System(vcat(waveguide_domains,sys.domains))
end
function add_pc_waveguide(sim::Simulation; x0=0, y0=0, direction, width)
    sys = add_defect_waveguide(sim.sys; x0=x0, y0=y0, direction=direction, width=width)
    return Simulation(sim; :sys => sys)
end
