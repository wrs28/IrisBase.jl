"""
    domain = Domain(; is_in_domain = whole_domain,
        domain_params = Dict(:n₁ => 1, :n₂ => 0, :F => 0),
        domain_type = :background,
        domain_ε::Function = piecewise_constant_ε,
        domain_F::Function = piecewise_constant_F,
        is_in_subdomain = Function[],
        subdomain_params = Dict{Symbol,Float64}[],
        subdomain_type = Symbol[],
        subdomain_ε = Function[],
        subdomain_F = Function[],
        lattice = Bravais(),
        which_asymptote = :none,
        which_waveguide = 0)

- `is_in_domain` is boolean function with inputs `x`, `y`, `domain_params`

INCOMPLETE DOCUMENTATION

domain = Domain(domain; :key1 => value1, :key2 => value2, ...)

New domain structure from old with modified fields
"""
struct Domain{TSH1,TDP,TDT,TDE,TDF,TSH2,TSP,TSD,TSE,TSF}
    is_in_domain::TSH1
    domain_params::TDP
    domain_type::TDT
    domain_ε::TDE
    domain_F::TDF

    is_in_subdomain::TSH2
    subdomain_params::TSP
    subdomain_type::TSD
    subdomain_ε::TSE
    subdomain_F::TSF

    lattice::BravaisLattice
    which_asymptote::Symbol
    which_waveguide::Int

    num_subdomains::Int

    function Domain(;
        is_in_domain::Tsh1 = Universe(),
        domain_params::Tdp = Dict(:n₁ => 1, :n₂ => 0, :F => 0),
        domain_type::Tdt = GenericDomain(),
        domain_ε::Tde = piecewise_constant_ε,
        domain_F::Tdf = piecewise_constant_F,
        is_in_subdomain::Tsh2 = (),
        subdomain_params::Tsp = (),
        subdomain_type::Tsd = (),
        subdomain_ε::Tse = (),
        subdomain_F::Tsf = (),
        lattice::BravaisLattice = BravaisLattice(),
        which_asymptote::Symbol = :none,
        which_waveguide::Int = 0
        ) where {Tsh1,Tdp,Tsh2<:Tuple,Tsp,Tsd,Tdt<:AbstractDomain,Tde<:Function,Tdf<:Function,Tse,Tsf}

        num_subdomains= 1 + length(is_in_subdomain)

        @assert !(which_asymptote==:none && which_waveguide!==0) "assigned waveguide id $(which_waveguide) to non-asymptotic domain"

        return new{Tsh1,Tdp,Tdt,Tde,Tdf,Tsh2,Tsp,Tsd,Tse,Tsf}(
                    is_in_domain, domain_params, domain_type, domain_ε, domain_F,
                    is_in_subdomain, subdomain_params, subdomain_type, subdomain_ε, subdomain_F,
                    lattice, which_asymptote, which_waveguide, num_subdomains
                )
    end

    function Domain(dom::Td;
        is_in_domain = dom.is_in_domain,
        domain_params = dom.domain_params,
        domain_type = dom.domain_type,
        domain_ε = dom.domain_ε,
        domain_F = dom.domain_F,
        is_in_subdomain = dom.is_in_subdomain,
        subdomain_params = dom.subdomain_params,
        subdomain_type = dom.subdomain_type,
        subdomain_ε = dom.subdomain_ε,
        subdomain_F = dom.subdomain_F,
        lattice = dom.lattice,
        which_asymptote = dom.which_asymptote,
        which_waveguide = dom.which_waveguide
        ) where Td<:Domain

        return Domain(;is_in_domain = is_in_domain, domain_params = domain_params, domain_type = domain_type,
                domain_ε = domain_ε, domain_F = domain_F, is_in_subdomain = is_in_subdomain, subdomain_params = subdomain_params,
                subdomain_type = subdomain_type, subdomain_ε = subdomain_ε, subdomain_F = subdomain_F, lattice = lattice,
                which_asymptote = which_asymptote, which_waveguide = which_waveguide)
    end

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
end


"""
    sys = System(domains::Array{Domain})
    sys = System(domain)

collection of domains that defines System, an input for Simulation.

input can be array or scalar.

    sys = System(sys)

system object from system object
"""
struct System{TDOM,TDIS,TPAR,TE,TF}
    domains::TDOM
    ε::Array{ComplexF64,2}
    F::Array{Float64,2}
    χ::TDIS

    domain_by_region::Array{Int,1}
    params_by_region::TPAR
    ε_by_region::TE
    F_by_region::TF

    regions::Array{Int,2}
    num_prev_regions::Array{Int,1}

    waveguides::Array{Int,1}

    function System(
                domains_in::Td = (Domain(),),
                χ::Tdis = TwoLevelSystem(),
                ε::Array{ComplexF64,2} = Array{ComplexF64}(undef, 0, 0),
                F::Array{Float64,2} = Array{Float64}(undef, 0, 0),
                regions::Array{Int,2} = Array{Int}(undef, 0, 0)
                ) where {Td<:Tuple,Tdis<:AbstractDispersion}

        domains1 = fsort(domains_in; by=isBackground)
        domains2 = bsort(domains1; by=isDefect)
        domains3 = fsort(domains2; by=isPC)
        domains = bsort(domains3; by=isWaveguide)

        num_regions = sum(map(x->x.num_subdomains,domains))

        params_by_region = join_tuple(map(d->(d.subdomain_params...,d.domain_params),domains)...)
        ε_by_region = join_tuple(map(d->(d.subdomain_ε...,d.domain_ε),domains)...)
        F_by_region = join_tuple(map(d->(d.subdomain_F...,d.domain_F),domains)...)
        domain_by_region = Array{Int}(undef, num_regions)
        num_prev_regions = collect(cumsum_tuple(map(x->x.num_subdomains,domains)))
        idx = 1
        for i ∈ eachindex(domains)
            for j ∈ 1:domains[i].num_subdomains-1
                domain_by_region[idx] = i ; idx += 1
            end
            domain_by_region[idx] = i ; idx += 1
        end
        waveguides = sort(unique(map(d->d.which_waveguide,domains)))[2:end]
        return new{Td,Tdis,typeof(params_by_region),typeof(ε_by_region),typeof(F_by_region)}(domains, ε, F, χ, domain_by_region, params_by_region, ε_by_region, F_by_region, regions, num_prev_regions, waveguides)
    end

    System(sys::Ts) where Ts<:System = System(sys.domains)

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
    end
end

join_tuple(x::Tuple, y...) = (x..., join_tuple(y...)...)
join_tuple(x::Tuple) = x

fsort(a::Tuple;kwargs...) =
    ssort(ssort(a...;kwargs...)...;kwargs...)
bsort(a::Tuple;kwargs...) =
    reverse(fsort(reverse(fsort(a;kwargs...));kwargs...))
ssort(a,b...;by) = by(a.domain_type) ? (a,ssort(b...;by=by)...) : (ssort(b...;by=by)...,a)
ssort(a;kwargs...)=(a,)

cumsum_tuple(a::Tuple) = reverse(_cumsum_tuple(reverse(a)...))
_cumsum_tuple(a,b...) = (sum(b),_cumsum_tuple(b...)...)
_cumsum_tuple(a) = 0


"""
    dis = Discretization(dis; key1 => value1, key2 => value2...)
new discretization object from old, with modified fields.
"""
struct Discretization{TCS}
    dx::NTuple{2,Float64}
    coordinate_system::TCS
    sub_pixel_num::Int
    origin::NTuple{2,Float64}

    N_tr::NTuple{2,Int}
    N::NTuple{2,Int}
    dN::NTuple{2,NTuple{2,Int}}

    x::NTuple{2,Array{Float64,2}}
    x_tr::NTuple{2,Array{Float64,2}}
    x_idx::NTuple{2,Array{Int,1}}

    X::Array{Float64,2}
    Y::Array{Float64,2}

    X_idx::Array{Int,1}

    function Discretization(dx::Tuple, coordinate_system::CS, sub_pixel_num::Int, origin::Tuple, N::Tuple, dN::NTuple{2,Tuple}) where CS

        N_tr = (N[1] - dN[1][1] - dN[1][2], N[2] - dN[2][1] - dN[2][2])

        x_tr = (Array{Float64}(undef,N_tr[1],1), Array{Float64}(undef,1,N_tr[2]))
        x = (Array{Float64}(undef,N[1],1), Array{Float64}(undef,1,N[2]))
        x_idx = (Array{Int}(undef, N_tr[1]), Array{Int}(undef, N_tr[2]))
        for j ∈ eachindex(N)
            if !isinf(dx[j])
                for i ∈ eachindex(x[j])
                    x[j][i] = origin[j] + dx[j]*(i+1/2)#(1/2:N[j]-1/2)
                end
                for i ∈ eachindex(x_tr[j])
                    x_tr[j][i] = x[j][1] + dN[j][1]*dx[j] + dx[j]*(i)#collect(0:N_tr[j]-1).*dx[j]
                end
                for i ∈ eachindex(x_idx[j])
                    x_idx[j][i] = dN[j][1] + i#(1:N_tr[j])
                end
            else
                x[j][:] .= origin[j]
                x_tr[j][:] .= x[j][1]
                x_idx[j][:] .= dN[j][1] .+ (1:N_tr[j])
            end
        end

        @assert CS ∈ (Cartesian,Polar) "unrecognized coordinate system type $CS"
        if CS == Cartesian
            X, Y = broadcast((x,y)->x,x[1],x[2]), broadcast((x,y)->y,x[1],x[2])
        else
            X = x[1].*cos.(x[2])
            Y = x[1].*sin.(x[2])
        end

        X_idx = LinearIndices(Array{Bool}(undef, N...))[x_idx...][:]
        return new{CS}(float.(dx), coordinate_system, sub_pixel_num, float.(origin), N_tr, N, dN,
            x, x_tr, x_idx, X, Y, X_idx)
    end

    Discretization(dis::Discretization{CS}; dx=dis.dx, coordinate_system=CS(), sub_pixel_num=dis.sub_pixel_num) where CS = Discretization(dx, coordinate_system, sub_pixel_num, (0.,0.), (1,1), ((0,0),(0,0)))

    function Base.show(io::IO, dis::Discretization)
        if !get(io, :sub, false)
            print(io, "Discretization: \n")
        end
        print(io, "\tN: ", dis.N, "\n",
        "\tsub-pixel number: ", dis.sub_pixel_num, "\n",
        "\tdx: (",fmt("1.4f",dis.dx[1]),",",fmt("1.4f",dis.dx[2]),")\n")
        if isCartesian(dis.coordinate_system)
            print(io, "\tcoordinates: Cartesian")
        elseif isPolar(dis.coordinate_system)
            print(io, "\tcoordinates: Polar")
        else
            print(io, "\tcoordinates: ", dis.coordinate_system)
        end
    end
end


"""
    bnd = Boundary(bnd; :key1 => value1, :key2 => value2, ...)

new boundary object from old, with modified fields.
"""
struct Boundary{TBC11,TBC12,TBC21,TBC22,TBL11,TBL12,TBL21,TBL22}
    ∂Ω::NTuple{2,NTuple{2,Float64}}
    ∂Ω_tr::NTuple{2,NTuple{2,Float64}}
    bc::Tuple{Tuple{TBC11,TBC12},Tuple{TBC21,TBC22}}
    bl::Tuple{Tuple{TBL11,TBL12},Tuple{TBL21,TBL22}}

    function Boundary(
                ∂Ω::NTuple{2,NTuple{2,Number}},
                BC::Tuple{Tuple{TBC1,TBC2},Tuple{TBC3,TBC4}} where {TBC1<:AbstractBC,TBC2<:AbstractBC,TBC3<:AbstractBC,TBC4<:AbstractBC},
                BL::Tuple{Tuple{TBL1,TBL2},Tuple{TBL3,TBL4}} where {TBL1<:AbstractBL,TBL2<:AbstractBL,TBL3<:AbstractBL,TBL4<:AbstractBL};
                warnings::Bool=true
                )
        if warnings
            (typeof(BL[1][1])<:Union{PML,cPML} && !(typeof(BC[1][1])<:DirichletBC)) ? (@warn "absorbing layer with non-dirichlet condition for dim 1, side 1, which is inconsistent") : nothing
            (typeof(BL[1][2])<:Union{PML,cPML} && !(typeof(BC[1][2])<:DirichletBC)) ? (@warn "absorbing layer with non-dirichlet condition for dim 1, side 2, which is inconsistent") : nothing
            (typeof(BL[2][1])<:Union{PML,cPML} && !(typeof(BC[2][1])<:DirichletBC)) ? (@warn "absorbing layer with non-dirichlet condition for dim 2, side 1, which is inconsistent") : nothing
            (typeof(BL[2][2])<:Union{PML,cPML} && !(typeof(BC[2][2])<:DirichletBC)) ? (@warn "absorbing layer with non-dirichlet condition for dim 2, side 2, which is inconsistent") : nothing
        end

        ∂Ω_tr1 = (∂Ω[1][1] + BL[1][1].depth, ∂Ω[1][2] - BL[1][2].depth)
        ∂Ω_tr2 = (∂Ω[2][1] + BL[2][1].depth, ∂Ω[2][2] - BL[2][2].depth)
        ∂Ω_tr = (∂Ω_tr1, ∂Ω_tr2)

        new{typeof(BC[1][1]),typeof(BC[1][2]),typeof(BC[2][1]),typeof(BC[2][2]),typeof(BL[1][1]),typeof(BL[1][2]),typeof(BL[2][1]),typeof(BL[2][2])}(map(x->map(float,x),∂Ω),map(x->map(float,x),∂Ω_tr),BC,BL)
    end

    Boundary(bnd::Boundary; ∂Ω=bnd.∂Ω, bc=bnd.bc, bl=bnd.bl) = Boundary(∂Ω, bc, bl)

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
            if T<:FloquetBC
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
end


"""
    channel = Channels(waveguide=0, quantum_number=0)

channel object for Scattering().

`waveguide` identifies which waveguide the channel is defined by

`quantum_number` identifies which mode of the waveguide defines the channel


    channel = Channels(chn; key1 => value1, key2 => value2, ...)

new channel object from old with modified fields.
"""
struct Channels
    waveguide::Int
    quantum_number::Int
    dispersion::Array{AbstractInterpolation,1}
    gaps::Array{Array{Float64,1},1}

    Channels(waveguide=0, quantum_number=0) = new(waveguide, quantum_number, AbstractInterpolation[], Array{Float64,1}[])

    Channels(chn::Channels; waveguide=chn.waveguide, quantum_number=chn.quantum_number) = Channels(waveguide, quantum_number)

    function Base.show(io::IO, chn::Channels)
        if get(io, :indented, false)
            print(io, "\t\twaveguide: ", chn.waveguide, "\n",
            "\t\tquantum number: ", chn.quantum_number)
        else
            print(io, "Channel: \n")
            print(io, "\twaveguide: ", chn.waveguide, "\n",
            "\tquantum number: ", chn.quantum_number, "\n")
        end
    end
end


"""
    sct = Scattering(channels::Array{Channel})

scattering object for Simulation

`channels` array of Channel objects


    sct = Scattering(sct; :key1 => value1, :key2 => value2, ...)

new scattering object from old with modified fields.
"""
struct Scattering
    channels::Array{Channels,1}
    waveguides_used::Array{Int,1}

    ε₀::Array{Array{ComplexF64,2},1}

    function Scattering(channels::Array{Channels,1}=Channels[], ε₀=[])

        channels = deepcopy(channels)
        ε₀ = deepcopy(ε₀)

        num_channels = length(channels)
        waveguides = Array{Int}(undef,num_channels)
        for i ∈ eachindex(waveguides)
            waveguides[i]=channels[i].waveguide
        end
        waveguides_used = unique(waveguides)
        ε₀ = Array{Array{ComplexF64,2},1}(undef,length(waveguides_used))
        return new(channels, waveguides_used, ε₀)
    end

    Scattering(sct::Scattering; channels=sct.channels) = Scattering(channels)

    function Base.show(io::IO, sct::Scattering)
        plot_flag=false
        if !get(io, :sub, false)
            plot_flag = true
            print(io, "Scattering with ", length(sct.channels), " channels:\n")
        end
        temp = [["\n", sct.channels[i]] for i ∈ eachindex(sct.channels)]
        for i ∈ eachindex(temp)
            if i>1
                print(io, "\n")
            end
            print(io, "\tChannel ", i, ":")
            print(IOContext(io::IO, :indented => true), temp[i]...)
        end
    end
end


"""
See also: [`System`](@ref), [`Boundary`](@ref), [`Discretization`](@ref), [`Scattering`](@ref), [`TwoLevelSystem`](@ref), [`Bravais`](@ref)

    sim = Simulation(; sys, bnd, dis, sct=Scattering(), tls=TwoLevelSystem(), lat=Bravais(bnd))

simulation object
"""
struct Simulation{TSYS,TBND,TDIS}
    sys::TSYS
    bnd::TBND
    dis::TDIS
    sct::Scattering
    lat::BravaisLattice

    function Simulation(
        bnd::Boundary,
        dis::Discretization{CS},
        sys::System,
        sct::Scattering,
        lat::BravaisLattice;
        display::Bool=false,
        warnings::Bool=true,
        smoothing::Bool=true
        ) where CS<:CoordinateSystem

        ∂Ω = bnd.∂Ω
        L = (∂Ω[1][2]-∂Ω[1][1],∂Ω[2][2]-∂Ω[2][1])

        dx = dis.dx
        @assert all(L./dx .< Inf) "lattice spacings $dx inconsistent with size $L"
        N = round.(Int,L./dx)

        dN = Array{Int}(undef,2,2)
        index1 = findmin(N)[2]
        index2 = mod1(index1+1,2)
        @assert !iszero(L[index2]) "no width along dimension $index2"
        if N[index1]≤1
            # N[index1] = 1
            # DX = [Inf, L[index2]/N[index2]]
            # dx = DX[[index1,index2]]
            # dN[:,index1] .= 0
            # dN[:,index2] .= floor.(Int, bnd.bl_depth[:,index2]/dx[index2])
            # for side ∈ 1:2
            #     bnd.BL[index2][side].depth = dN[side,index2]*dx[index2]
            # end
        else
            DX = [L[index1]/N[index1], L[index2]/round(Int,L[index2]/dx[index2])]
            dx = (DX[index1],DX[index2])
            N = round.(Int,L./dx)
            dN = ()
            for dim ∈ 1:2
                dNside = ()
                for side ∈ 1:2
                    dNside = (dNside...,floor(Int, bnd.bl[dim][side].depth/dx[dim]))
                end
                dN = (dN...,dNside)
            end

            for dim ∈ 1:2, side ∈ 1:2
                bnd.bl[dim][side].depth = dN[dim][side]*dx[dim]
            end
        end

        @assert all(sum.(dN) .< N) "boundary layers fill entire region. increase size of domain or decrease size of boundary layers."

        origin = (bnd.∂Ω[1][1],bnd.∂Ω[2][1])
        if typeof(bnd.bc[1][1])<:MatchedBC || typeof(bnd.bc[1][2])<:MatchedBC
            if typeof(bnd.bc[2][1])<:MatchedBC || typeof(bnd.bc[2][2])<:MatchedBC
                throw()
            else
                non_matched_dim = 2
            end
        else
            if typeof(bnd.bc[2][1])<:MatchedBC || typeof(bnd.bc[2][2])<:MatchedBC
                throw()
            else
                non_matched_dim = 1
            end
        end

        bb = apply_args(bnd.bc, CS(), bnd.bc[non_matched_dim][1], bnd.bc[non_matched_dim][2]; N=N, dx=dx, xmin=origin, lattice=lat)
        bnd = Boundary(∂Ω, bb, apply_args(bnd.bl; Δ=∂Ω))
        dis = Discretization(dx, dis.coordinate_system, dis.sub_pixel_num, origin, N, dN)

        # argument checking
        for side ∈ 1:2
            if warnings
                for dim ∈ 1:2
                    typeof(bnd.bc[dim][side])<:FloquetBC ? (typeof(bnd.bc[dim][mod1(side+1,2)])<:FloquetBC ? nothing : (@warn "only one side=$side of dim=$dim set to floquet")) : nothing
                end
            end
            @assert ((typeof(bnd.bc[1][side])<:FloquetBC && !isinf(lat.a)) || !(typeof(bnd.bc[1][side])<:FloquetBC)) "floquet bc in dim=1, but infinite lattice constant a"
            @assert ((typeof(bnd.bc[2][side])<:FloquetBC && !isinf(lat.b)) || !(typeof(bnd.bc[1][side])<:FloquetBC)) "floquet bc in dim=2, but infinite lattice constant b"
        end
        if warnings
            for i ∈ eachindex(sys.domains)
                if N[1] == 1 && !isinf(sys.domains[i].lattice.a)
                    @warn "one-dimensional (vertical) system, but domain $(i) has finite transverse lattice constant $(sys.domains[i].lattice.a) < ∞"
                end
                if N[2] == 1 && !isinf(sys.domains[i].lattice.b)
                    @warn "one-dimensional (horizontal) system, but domain $(i) has finite transverse lattice constant $(sys.domains[i].lattice.b) < ∞"
                end
            end
        end

        ε, F, regions = construct_εFr(bnd,dis,sys)
        smoothing ? sub_pixel_smoothing!(bnd, dis, sys, ε, F, regions; display=display) : nothing

        sys = System(sys.domains, sys.χ, ε, F, regions)

        which_waveguides = Array{Int,1}(undef,length(sys.domains))
        for j ∈ eachindex(sys.domains)
            which_waveguides[j] = sys.domains[j].which_waveguide
        end
        for i ∈ 1:length(sct.waveguides_used)
            sct.ɛ₀[i] = Array{ComplexF64}(undef,dis.N[1],dis.N[2])
            wg_inds = findall(which_waveguides .== sct.waveguides_used[i])
            wg_domains = Array{Domain}(undef,length(wg_inds))
            for j ∈ eachindex(wg_domains)
                wg_domains[j] = Domain(sys.domains[wg_inds[j]]; :which_asymptote => :none, :which_waveguide => 0)
            end
            if isempty(wg_domains)
                @warn "channels references undefined waveguide $(sct.waveguides_used[i]), scattering calculations will fail"
            else
                sct.ɛ₀[i][:] = sub_pixel_smoothing(bnd, dis, System(wg_domains), display=display)[1][:]
            end
        end
        return new{typeof(sys), typeof(bnd), typeof(dis)}(sys, bnd, dis, sct, lat)
    end

    Simulation(sim::Simulation; sys=sim.sys, bnd=sim.bnd, dis=sim.dis, sct=sim.sct, tls=sim.tls, lat=sim.lat, disp_opt=false) = deepcopy(Simulation(sys=System(sys), bnd=Boundary(bnd), dis=Discretization(dis), sct=Scattering(sct), tls=TwoLevelSystem(tls), lat=BravaisLattice(lat), disp_opt=disp_opt))

    function Base.show(io::IO, sim::Simulation)
        print(IOContext(io, :sub=>true),
        "Simulation : \n\n",
        "sys: \n", sim.sys, "\n\n",
        "bnd: \n", sim.bnd, "\n\n",
        "dis: \n", sim.dis, "\n\n",
        "sct: \n", sim.sct)
        if [sim.lat.a, sim.lat.b] !== [Inf,Inf]
            print(IOContext(io, :sub=>true, :sub1=>false), "\n\nlat: ", sim.lat)
        end
    end
end
