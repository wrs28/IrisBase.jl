#
# is_in_subdomain::Array{TS2,1}
# subdomain_params::Array{TSP,1}
# subdomain_type::Array{Symbol,1}
# subdomain_ε::Array{Function,1}
# subdomain_F::Array{Function,1}
# num_subdomains::Int

### DOMAIN STRUCT
function Base.show(io::IO, dom::Domain)
    print(io, "Domain: \n")
    print(io, "\tDomain type: ", dom.domain_type, "\n")
    if dom.which_asymptote !== :none
        print(io, "\t\twaveguide: ", dom.which_waveguide,"\n")
        print(io, "\t\tasymptotic region: ", dom.which_asymptote, "\n")
    end
    print(io, "\tdomain shape: ", dom.is_in_domain, "\n")
    print(IOContext(io, :typeinfo => Dict), "\tdomain params: ", dom.domain_params, "\n")
    print(IOContext(io, :sub=>true, :sub1=>true), "\tdomain lattice: ", dom.lattice, "\n")
    print(IOContext(io, :typeinfo => Array{Function}), "\tbackground dielectric function: ", dom.domain_ε, "\n")
    print(IOContext(io, :typeinfo => Array{Float64}), "\tbackground pump function : ", dom.domain_F, "\n")
        print(IOContext(io, :typeinfo => Array{AbstractShape}), "\tsubdomain shapes: ", dom.is_in_subdomain, "\n")
    print(IOContext(io, :typeinfo => Array{Dict{Symbol},1}), "\tsubdomain params: ", dom.subdomain_params)
end


### SYSTEM STRUCT
function Base.show(io::IO, sys::System)
    if !get(io, :sub, false)
        print(io, "System with ", length(sys.domains), " domains: \n")
    end
    domain_string = [["\tdomain ", i, " type: ", sys.domains[i].domain_type] for i ∈ eachindex(sys.domains)]
    asymptote_string = [[", ", sys.domains[i].which_asymptote] for i ∈ eachindex(sys.domains)]
    waveguide_string = [[", waveguide ", sys.domains[i].which_waveguide, "\n"] for i ∈ eachindex(sys.domains)]
    for i ∈ eachindex(sys.domains)
        print(io, domain_string[i]...)
        if sys.domains[i].which_asymptote !== :none
            print(io, asymptote_string[i]...)
            print(io, waveguide_string[i]...)
        else
            print(io, "\n")
        end
    end
    for i ∈ eachindex(sys.waveguides)
        if i==1
            print(io,"\n")
        end
        print(IOContext(io, :typeinfo => Array{Int}), "\twaveguide ", sys.waveguides[i], " domains: ", findall(sys.waveguides[i].==[sys.domains[j].which_waveguide for j ∈ eachindex(sys.domains)]), "\n")
        if i==length(sys.waveguides)
            print(io,"\n")
        end
    end
    # print(IOContext(io, :typeinfo => Array{Int}), "\tdomain: ", sys.domain_by_region, "\n")
    # print(IOContext(io, :typeinfo => Array{ComplexF64}), "\tn: ", sys.n_by_region, "\n")
    # print(IOContext(io, :typeinfo => Array{ComplexF64}), "\tε: ", sys.ε_by_region, "\n")
    # print(IOContext(io, :typeinfo => Array{Float64}), "\tF: ", sys.F_by_region)
end


### BOUNDARY STRUCT
Base.show(io::IO, bnd::Boundary) = begin
    if !get(io, :sub, false)
        print(io, "Boundary: \n")
    end

    print(io, "\t\t\t-------------------------\n")

    print(io, "\tbound.\t|\t\t   ", fmt("6s",_get_b_symbol(bnd.bc[2][2])), "\t\t|\n")
    print(io, "\tcond-\t|   ",fmt("6s",_get_b_symbol(bnd.bc[1][1])), "\t\t   ", fmt("5s",_get_b_symbol(bnd.bc[1][2])), "|\n")
    print(io, "\t ition\t|\t\t   ", fmt("6s",_get_b_symbol(bnd.bc[2][1])), "\t\t|\n")

    print(io, "\t\t\t-------------------------\n")

    print(io, "\t\t\t|\t\t", fmt("+5.3f",bnd.∂Ω[2][2]), "\t\t\t|\n")
    print(io, "\t∂Ω\t\t| ", fmt("+5.3f",bnd.∂Ω[1][1]), "\t\t", fmt("+5.3f",bnd.∂Ω[1][2]), "\t|\n")
    print(io, "\t\t\t|\t\t", fmt("+5.3f",bnd.∂Ω[2][1]), "\t\t\t|\n")

    print(io, "\t\t\t-------------------------\n")

    print(io, "\t\t\t|\t\t", fmt("+5.3f",bnd.∂Ω_tr[2][2]), "\t\t\t|\n")
    print(io, "\t∂Ω_tr\t| ",fmt("+5.3f",bnd.∂Ω_tr[1][1]), "\t\t", fmt("+5.3f",bnd.∂Ω_tr[1][2]), "\t|\n")
    print(io, "\t\t\t|\t\t", fmt("+5.3f",bnd.∂Ω_tr[2][1]), "\t\t\t|\n")

    print(io, "\t\t\t-------------------------\n")

    print(io, "\tbound.\t|\t\t", fmt("8s",_get_b_symbol(bnd.bl[2][2])), "\t\t|\n")
    print(io, "\tlayer\t|",fmt("8s",_get_b_symbol(bnd.bl[1][1])), "\t\t", fmt("8s",_get_b_symbol(bnd.bl[1][2])), "|\n")
    print(io, "\t\t\t|\t\t", fmt("8s",_get_b_symbol(bnd.bl[2][1])), "\t\t|\n")

    print(io, "\t\t\t-------------------------\n")

    print(io, "\tbound.\t|\t\t", fmt("+5.3f",bnd.bl[2][2].depth), "\t\t\t|\n")
    print(io, "\tlayer\t| ",fmt("+5.3f",bnd.bl[1][1].depth), "\t\t", fmt("+5.3f",bnd.bl[1][2].depth), "\t|\n")
    print(io, "\tdepth\t|\t\t", fmt("+5.3f",bnd.bl[2][1].depth), "\t\t\t|\n")

    print(io, "\t\t\t-------------------------\n")
end
function _get_b_symbol(bcl::T) where T<:Union{AbstractBC,AbstractBL}
    if T<:AbstractBC
        if T<:PeriodicBC
            return :p
        elseif T<:RobinBC
            return :r
        elseif T<:DirichletBC
            return :d
        elseif T<:NeumannBC
            return :n
        elseif T<:MatchedBC
            return :m
        else
            throw(ArgumentError("unrecognized boundary condition"))
        end
    else
        if T<:PML
            return :PML
        elseif T<:cPML
            return :cPML
        elseif T<:noBL
            return :none
        else
            throw(ArgumentError("unrecognized boundary layer"))
        end
    end
end


### DISCRETIZATION STRUCT
Base.show(io::IO, dis::Discretization) = begin
    if !get(io, :sub, false)
        print(io, "Discretization: \n")
    end
    print(io, "\tN: ", dis.N, "\n",
    # "\ttruncated N: ", dis.N_tr, "\n",
    "\tsub-pixel number: ", dis.sub_pixel_num, "\n",
    "\tdx: ", dis.dx, "\n")
    if isCartesian(dis.coordinate_system)
        print(io, "\tcoordinates: Cartesian")
    elseif isPolar(dis.coordinate_system)
        print(io, "\tcoordinates: Polar")
    else
        print(io, "\tcoordinates: ", dis.coordinate_system)
    end
end


### CHANNEL STRUCT
Base.show(io::IO, chn::Channels) = begin
    if get(io, :indented, false)
        print(io, "\t\twaveguide: ", chn.waveguide, "\n",
        "\t\tquantum number: ", chn.quantum_number)
    else
        print(io, "Channel: \n")
        print(io, "\twaveguide: ", chn.waveguide, "\n",
        "\tquantum number: ", chn.quantum_number, "\n")
        # if !isempty(chn.dispersion)
        #     dispersion_plot = lineplot(x->chn.dispersion[1](x), 0:.01:2π, canvas=BlockCanvas, color=:yellow, title="Dispersion", name="channel");
        #     for j ∈ eachindex(chn.gaps)
        #         lineplot!(dispersion_plot, chn.gaps[1][1], 0, color=:white);
        #         lineplot!(dispersion_plot, chn.gaps[1][2], 0, color=:white);
        #     end
        #     println(dispersion_plot)
        # end
    end
end


### SCATTERING STRUCT
Base.show(io::IO, sct::Scattering) = begin
    plot_flag=false
    if !get(io, :sub, false)
        plot_flag = true
        print(io, "Scattering with ", length(sct.channels), " channels:\n")
    end
    temp = [["\n", sct.channels[i]] for i ∈ eachindex(sct.channels)]
    # if plot_flag && !isempty(sct.channels) && !isempty(sct.channels[1].dispersion)
        # dispersion_plot = lineplot(x->sct.channels[1].dispersion[1](x), sct.channels[1].dispersion[1].ranges[1][1], sct.channels[1].dispersion[1].ranges[1][end], canvas=BlockCanvas, title="Dispersion", name="channel 1")
    # end
    for i ∈ eachindex(temp)
        if i>1
            print(io, "\n")
        end
        print(io, "\tChannel ", i, ":")
        print(IOContext(io::IO, :indented => true), temp[i]...)
        # if plot_flag
            # if !isempty(sct.channels[i].dispersion)
                # if i>1
                    # lineplot!(dispersion_plot,x->sct.channels[i].dispersion[1](x), sct.channels[i].dispersion[1].ranges[1][1], sct.channels[i].dispersion[1].ranges[1][end], name="channel $i");
                # end
                # for j ∈ eachindex(sct.channels[i].gaps)
                    # lineplot!(dispersion_plot, [0,1], fill(sct.channels[i].gaps[1][1],2), color=:white);
                    # lineplot!(dispersion_plot, sct.channels[i].gaps[1][2], 0, color=:white);
                # end
            # end
        # end
    end
    # try
        # dispersion_plot
        # print(io,"\n\n")
        # println(dispersion_plot)
    # catch
        # nothing
    # end
end


### TWO-LEVEL-SYSTEM STRUCT
Base.show(io::IO, tls::TwoLevelSystem) = begin
    if !get(io, :sub, false)
        print(io, "Two Level System:\n")
    end
    print(io, "\tD₀: ", tls.D₀, "\n",
    "\tω₀: ", tls.ω₀, "\n",
    "\tγ⟂: ", tls.γp)
end


### SIMULATION STRUCT
Base.show(io::IO, sim::Simulation) = begin
    print(IOContext(io, :sub=>true),
    "Simulation : \n\n",
    "sys: \n", sim.sys, "\n\n",
    "bnd: \n", sim.bnd, "\n\n",
    "dis: \n", sim.dis, "\n\n",
    "sct: \n", sim.sct, "\n\n",
    "tls: \n", sim.tls)
    if [sim.lat.a, sim.lat.b] !== [Inf,Inf]
        print(IOContext(io, :sub=>true, :sub1=>false), "\n\nlat: ", sim.lat)
    end
end
