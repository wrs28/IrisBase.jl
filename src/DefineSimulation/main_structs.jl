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
"""
struct Domain{TS1,TDP,TS2,TSP}
    is_in_domain::TS1
    domain_params::TDP
    domain_type::Symbol
    domain_ε::Function
    domain_F::Function

    is_in_subdomain::Array{TS2,1}
    subdomain_params::Array{TSP,1}
    subdomain_type::Array{Symbol,1}
    subdomain_ε::Array{Function,1}
    subdomain_F::Array{Function,1}

    lattice::BravaisLattice
    which_asymptote::Symbol
    which_waveguide::Int

    num_subdomains::Int

    function Domain(;
        is_in_domain::Tsh1 = Universe(),
        domain_params::Tdp = Dict(:n₁ => 1, :n₂ => 0, :F => 0),
        domain_type::Symbol = :background,
        domain_ε::Function = piecewise_constant_ε,
        domain_F::Function = piecewise_constant_F,
        is_in_subdomain::Array{Tsh2,1} = AbstractShape[],
        subdomain_params::Array{Tsp,1} = Dict{Symbol,Float64}[],
        subdomain_type::Array{Symbol,1} = Symbol[],
        subdomain_ε::Array{Function} = Function[],
        subdomain_F::Array{Function} = Function[],
        lattice::BravaisLattice = BravaisLattice(),
        which_asymptote::Symbol = :none,
        which_waveguide::Int = 0
        ) where {Tsh1,Tdp,Tsh2,Tsp}

        num_subdomains= 1 + length(is_in_subdomain)

        @assert !(which_asymptote==:none && which_waveguide!==0) "assigned waveguide id $(which_waveguide) to non-asymptotic domain"

        return new{Tsh1,Tdp,Tsh2,Tsp}(
                    is_in_domain, domain_params, domain_type, domain_ε, domain_F,
                    is_in_subdomain, subdomain_params, subdomain_type, subdomain_ε, subdomain_F,
                    lattice, which_asymptote, which_waveguide, num_subdomains
                )
    end
end


"""
    sys = System(domains::Array{Domain})
    sys = System(domain)

collection of domains that defines System, an input for Simulation.

input can be array or scalar.
"""
struct System{TDOM,TPAR,TDIS}
    domains::Array{TDOM,1}
    ε::Array{ComplexF64,2}
    F::Array{Float64,2}
    χ::TDIS

    domain_by_region::Array{Int,1}
    params_by_region::Array{TPAR,1}
    ε_by_region::Array{Function,1}
    F_by_region::Array{Function,1}

    regions::Array{Int,2}
    num_prev_regions::Array{Int,1}

    waveguides::Array{Int,1}

    function System(
                domains::Array{Td,1} = [Domain()],
                χ::Tdis = TwoLevelSystem(),
                ε::Array{ComplexF64,2} = Array{ComplexF64}(undef, 0, 0),
                F::Array{Float64,2} = Array{Float64}(undef, 0, 0),
                regions::Array{Int,2} = Array{Int}(undef, 0, 0)
                ) where {Td<:Domain,Tdis<:AbstractDispersion}

        # order domains
        sort!(domains, by=isBulkWaveguide, rev=false)
        sort!(domains, by=isBackground, rev=false)
        sort!(domains, by=isDefect, rev=true)
        sort!(domains, by=isPC, rev=false)
        sort!(domains, by=isWaveguide, rev=true)

        num_regions = 0
        for d ∈ domains
            num_regions += d.num_subdomains
        end

        domain_by_region = Array{Int}(undef, num_regions)
        params_by_region = Array{typeof(domains[1].domain_params)}(undef, num_regions)
        ε_by_region = Array{Function}(undef, num_regions)
        F_by_region = Array{Function}(undef, num_regions)

        num_prev_regions = Array{Int}(undef, length(domains))
        num_prev_regions[1] = 0

        idx = 1
        for i ∈ eachindex(domains)
            d = domains[i]
            if i > 1
                num_prev_regions[i] = domains[i-1].num_subdomains
            end
            for j ∈ 1:d.num_subdomains-1
                domain_by_region[idx] = i
                params_by_region[idx] = d.subdomain_params[j]
                ε_by_region[idx] = d.subdomain_ε[j]
                F_by_region[idx] = d.subdomain_F[j]
                idx += 1
            end
            domain_by_region[idx] = i
            params_by_region[idx] = d.domain_params
            ε_by_region[idx] = d.domain_ε
            F_by_region[idx] = d.domain_F
            idx += 1
        end
        cumsum!(num_prev_regions, num_prev_regions, dims=1)

        waveguides = sort(unique([domains[i].which_waveguide for i ∈ eachindex(domains)]))[2:end]

        return new{Td,eltype(params_by_region),Tdis}(domains, ε, F, χ, domain_by_region, params_by_region, ε_by_region, F_by_region, regions, num_prev_regions, waveguides)
    end
end


struct Discretization{TCS}
    dx::NTuple{2,Float64}
    coordinate_system::Type{TCS}
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

    function Discretization(dx::Tuple, coordinate_system::Type, sub_pixel_num::Int, origin::Tuple, N::Tuple, dN::NTuple{2,Tuple})

        N_tr = (N[1] - dN[1][1] - dN[1][2], N[2] - dN[2][1] - dN[2][2])

        x_tr = (Array{Float64}(undef,N_tr[1],1), Array{Float64}(undef,1,N_tr[2]))
        x = (Array{Float64}(undef,N[1],1), Array{Float64}(undef,1,N[2]))
        x_idx = (Array{Int}(undef, N_tr[1]), Array{Int}(undef, N_tr[2]))
        for j ∈ eachindex(N)
            if !isinf(dx[j])
                x[j][:] .= origin[j] .+ dx[j]*(1/2:N[j]-1/2)
                x_tr[j][:] .= x[j][1] .+ dN[j][1]*dx[j] .+ collect(0:N_tr[j]-1).*dx[j]
                x_idx[j][:] .= dN[j][1] .+ (1:N_tr[j])
            else
                x[j][:] .= origin[j]
                x_tr[j][:] .= x[j][1]
                x_idx[j][:] .= dN[j][1] .+ (1:N_tr[j])
            end
        end

        @assert coordinate_system ∈ [Cartesian,Polar] "unrecognized coordinate system type $coordinate_system"
        if coordinate_system == Cartesian
            X, Y = broadcast((x,y)->x,x[1],x[2]), broadcast((x,y)->y,x[1],x[2])
        else
            X = x[1].*cos.(x[2])
            Y = x[1].*sin.(x[2])
        end

        X_idx = LinearIndices(Array{Bool}(undef, N...))[x_idx...][:]
        return new{coordinate_system}(float.(dx), coordinate_system, sub_pixel_num, float.(origin), N_tr, N, dN,
            x, x_tr, x_idx, X, Y, X_idx)
    end
end


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
end


"""
    channel = Channels(waveguide=0, quantum_number=0)

channel object for Scattering().

`waveguide` identifies which waveguide the channel is defined by

`quantum_number` identifies which mode of the waveguide defines the channel
"""
struct Channels
    waveguide::Int
    quantum_number::Int
    dispersion::Array{AbstractInterpolation,1}
    gaps::Array{Array{Float64,1},1}

    Channels(waveguide=0, quantum_number=0) = new(waveguide, quantum_number, AbstractInterpolation[], Array{Float64,1}[])
end


"""
    sct = Scattering(channels::Array{Channel})

scattering object for Simulation

`channels` array of Channel objects
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
    tls::TwoLevelSystem
    lat::BravaisLattice

    function Simulation(
        bnd::Boundary,
        dis::Discretization{CS},
        sys::System,
        sct::Scattering,
        tls::TwoLevelSystem,
        lat::BravaisLattice;
        display::Bool=false,
        warnings::Bool=true,
        smoothing::Bool=true
        ) where CS<:CoordinateSystem

        # shield structures from downstream changes
        # sys = deepcopy(sys); bnd = deepcopy(bnd); dis = deepcopy(dis)
        # sct = deepcopy(sct); tls = deepcopy(tls); lat = deepcopy(lat)

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
        bb = map(x->map(<:,x,(MatchedBC,MatchedBC)),map(x->map(typeof,x),bnd.bc))
        if true ∈ (bb[1]...,bb[2]...)
            if true ∈ bb[1] && true ∈ bb[2] && false ∉ (bb[1]...,bb[2])
            elseif true ∈ bb[1]
                non_matched_dim = 2
            elseif true ∈ bb[2]
                non_matched_dim = 1
            end
        else
            non_matched_dim = 0
        end

        bb = apply_args(bnd.bc, CS, typeof(bnd.bc[non_matched_dim][1]), typeof(bnd.bc[non_matched_dim][2]); N=N, dx=dx, xmin=origin, lattice=lat)
        bnd = Boundary(∂Ω, bb, apply_args(bnd.bl; Δ=∂Ω))
        dis = Discretization(dx, dis.coordinate_system, dis.sub_pixel_num, origin, N, dN)

        # argument checking
        for side ∈ 1:2
            if warnings
                for dim ∈ 1:2
                    typeof(bnd.bc[dim][side])<:PeriodicBC ? (typeof(bnd.bc[dim][mod1(side+1,2)])<:PeriodicBC ? nothing : (@warn "only one side=$side of dim=$dim set to periodic")) : nothing
                end
            end
            @assert ((typeof(bnd.bc[1][side])<:PeriodicBC && !isinf(lat.a)) || !(typeof(bnd.bc[1][side])<:PeriodicBC)) "periodic bc in dim=1, but infinite lattice constant a"
            @assert ((typeof(bnd.bc[2][side])<:PeriodicBC && !isinf(lat.b)) || !(typeof(bnd.bc[1][side])<:PeriodicBC)) "periodic bc in dim=2, but infinite lattice constant b"
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
        if smoothing
            ɛ, F = sub_pixel_smoothing(bnd, dis, sys, ε, F, regions; display=display)
        end

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

        return new{typeof(sys), typeof(bnd), typeof(dis)}(sys, bnd, dis, sct, tls, lat)
    end
end
