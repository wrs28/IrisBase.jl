################################################################################
### DOMAIN
################################################################################
"""
    domain = Domain(domain; :key1 => value1, :key2 => value2, ...)

New domain structure from old with modified fields
"""
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

################################################################################
### SYSTEM
################################################################################
"""
    sys = System(domain1, domain2, ...)
"""
function System(args::Vararg{Domain,N}) where N
    return System(vcat(args...))
end
"""
    sys = System(sys)

system object from system object
"""
function System(sys::Ts) where Ts<:System
    return System(sys.domains)
end

################################################################################
### DISCRETIZATION
################################################################################
# get_coordinate_system(dis::Discretization{CS}) where CS<:CoordinateSystem = CS

"""
    dis = Discretization{coordinate_system=Cartesian}(dx[, sub_pixel_num=DEFAULT_SUBSAMPLE_NUMBER])

discretization object for Simulation

`dx` lattice spacing (scalar)

`sub_pixel_num` is the subsampling rate used in sub-pixel smoothing.
"""
Discretization(dx::Union{Number,Array,Tuple}) = Discretization(dx, Cartesian)
Discretization(dx::Real, coordinate_system::Type{T}, sub_pixel_num::Int=DEFAULT_SUBSAMPLE_NUMBER) where T<:CoordinateSystem = Discretization((dx,dx), coordinate_system, sub_pixel_num, (0.,0.), (1,1), ((0,0),(0,0)))
Discretization(dx::Union{Array,Tuple}, coordinate_system::Type{T}, sub_pixel_num::Int=DEFAULT_SUBSAMPLE_NUMBER) where T<:CoordinateSystem = Discretization(Tuple(dx), coordinate_system, sub_pixel_num, (0.,0.), (1,1), ((0,0),(0,0)))
"""
    dis = Discretization(dis; key1 => value1, key2 => value2...)

new discretization object from old, with modified fields.
"""
Discretization(dis::Discretization{CS}; dx=dis.dx, coordinate_system=CS, sub_pixel_num=dis.sub_pixel_num) where CS = Discretization(dx, coordinate_system, sub_pixel_num, (0.,0.), (1,1), ((0,0),(0,0)))

################################################################################
### BOUNDARY
################################################################################
Boundary(∂Ω, bcls...) = begin; println(typeof.(_bnd_bcls(bcls))); end#Boundary(_bnd_Ω(∂Ω),_bnd_bcls(bcls)...); end

_bnd_Ω(∂Ω::NTuple{2,NTuple{2,Number}}) = ∂Ω
_bnd_Ω(∂Ω::NTuple{4,Number}) = ( (∂Ω[1],∂Ω[2]), (∂Ω[3],∂Ω[4]) )
_bnd_Ω(∂Ω::Array{T,1}) where T<:Number = ( (∂Ω[1],∂Ω[2]), (∂Ω[3],∂Ω[4]) )
_bnd_Ω(∂Ω::Array{Array{T,1},1}) where T<:Number = ( (∂Ω[1][1],∂Ω[1][2]), (∂Ω[2][1],∂Ω[2][2]) )
_bnd_Ω(∂Ω::Array{T,2}) where T<:Number = ( (∂Ω[1,1],∂Ω[1,2]), (∂Ω[2,1],∂Ω[2,2]) )

_bnd_bcls(bcs::NTuple{2,Tuple},bls::NTuple{2,Tuple}) = (bcs,bls)
_bnd_bcls(bcs::Tuple,bls::NTuple{2,Tuple}) = _bnd_bcls(bcs,bls...)
_bnd_bcls(bcs::NTuple{2,Tuple},bls::Tuple) = _bnd_bcls(bcs...,bls)
_bnd_bcls(x::Tuple,y::Tuple,bcls::Vararg{Tuple,N}) where N = _bnd_bcls((x...,y...),bcls...)
_bnd_bcls(bcls::Tuple) = _bnd_bcls(bcls...)
function _bnd_bcls(bcls::Vararg{Union{AbstractBC,AbstractBL},N}) where {N}
     bcs = bcls[map(<:, typeof.(bcls), fill(AbstractBC,N))]
     @assert length(bcs)==4 "insufficient number of boundaries ($(length(bcs))) specified"
     bcs = reorder(bcs)
println(bcs)
     bls = bcls[map(<:, typeof.(bcls), fill(AbstractBL,N))]
     bls = _oc_bls(bls)
     bls = reorder(bls)

     return ( (bcs[1],bcs[2]), (bcs[3],bcs[4]) ), ( (bls[1],bls[2]), (bls[3],bls[4]) )
end
"""
    bnd = Boundary(;∂Ω, bc=:d, bl=:none, bl_depth=0)

boundary object for Simulation.

`∂Ω` array of boundaries of computational domain. If not given, array of NaNs

`bc` scalar or array of boundary conditions

`bl` scalar or array of boundary layer types

`bl_depth` scalar or array of boundary layer depths
# """
function Boundary(;∂Ω=fill(NaN,4), bc::Type=DirichletBC, bl=noBL, depth::Real=0, r::Real=NaN, kwargs...)
    if !isnan(r)
        @assert bc<:AbstractBC "type of bc is $bc. when radius is specified, bc must be a boundary condition"
        @assert bl<:AbstractBL "type of bl is $bl. when radius is specified, bl must be a boundary layer"
        ∂Ω = (0, r, 0, 2π)
        if bc<:MatchedBC
            BC = bc{1,2}(kwargs[:outgoing_qns],kwargs[:incoming_qns])
        else
            BC = bc{1,2}(kwargs...)
        end
        BCs = ( DirichletBC{1,1}(), BC, PeriodicBC{2,1}(BravaisLattice(b=2π)), PeriodicBC{2,2}(BravaisLattice(b=2π)) )
        BLs = ( noBL{1,1}(0), bl{1,2}(depth), noBL{2,1}(0), noBL{2,2}(0) )
    else
        if bc == DirichletBC
            BCs = (DirichletBC{1,1}(),DirichletBC{1,2}(),DirichletBC{2,1}(),DirichletBC{2,2}())
        elseif bc == NeumannBC
            BCs = (NeumannBC{1,1}(),NeumannBC{1,2}(),NeumannBC{2,1}(),NeumannBC{2,2}())
        elseif bc == PeriodicBC
            BCs = (PeriodicBC{1,1}(),PeriodicBC{1,2}(),PeriodicBC{2,1}(),PeriodicBC{2,2}())
        else
            throw(ArgumentError("$(bc) incompatible with this wrapper signature of Boundary. For more control, use another signature"))
        end

        if bl == PML
            BLs = (PML{1,1}(depth),PML{1,2}(depth),PML{2,1}(depth),PML{2,2}(depth))
        elseif bl == cPML
            BLs = (cPML{1,1}(depth),cPML{1,2}(depth),cPML{2,1}(depth),cPML{2,2}(depth))
        elseif bl == noBL
            BLs = (noBL{1,1}(),noBL{1,2}(),noBL{2,1}(),noBL{2,2}())
        else
            throw(ArgumentError("unrecognized bl $(bl)"))
        end
    end

    return Boundary(∂Ω, BCs..., BLs...)
end
"""
    bnd = Boundary(bnd; :key1 => value1, :key2 => value2, ...)

new boundary object from old, with modified fields.
"""
Boundary(bnd::Boundary; ∂Ω=bnd.∂Ω, bc=bnd.bc, bl=bnd.bl) = Boundary(∂Ω, bc, bl)

################################################################################
### CHANNELS
################################################################################
"""
    channel = Channels(chn; key1 => value1, key2 => value2, ...)

new channel object from old with modified fields.
"""
Channels(chn::Channels; waveguide=chn.waveguide, quantum_number=chn.quantum_number) = Channels(waveguide, quantum_number)

################################################################################
### SCATTERING
################################################################################
Scattering(channel::Channels) = Scattering([channel])
"""
    sct = Scattering(sct; :key1 => value1, :key2 => value2, ...)

new scattering object from old with modified fields.
"""
Scattering(sct::Scattering; channels=sct.channels) = Scattering(channels)
"""
    sct = Scattering(channel1, channel2, ...)
"""
Scattering(args::Vararg{Channels}) = Scattering(vcat(args...))

################################################################################
### TWOLEVELSYSTEM
################################################################################
"""
    tls = TwoLevelSystem(tls; :key1 => value1, :key2 => value2, ...)

new tls object from old, with modified fields
"""
TwoLevelSystem(tls::TwoLevelSystem; D₀=tls.D₀, k₀=tls.k₀, γp=tls.γp, ω₀=tls.ω₀) = TwoLevelSystem(tls.D₀, tls.k₀, tls.γp, tls.ω₀)

################################################################################
### SIMULATION
################################################################################
"""
    sim = Simulation(sim; key1 => val1, key2 => val2, ...)

new simulation object from old, with modified fields.
"""
find_type(t::Type{T},x::T,y...) where T = x
find_type(t::Type{T},x, y...) where T = find_type(t,y...)
find_type(::Type{T},x) where T = nothing
function Simulation(args::Vararg{Any,N}; kwargs...) where N

    types = typeof.(args)
    simargs = ()

    bnd_ind = findfirst(map(<:,types,fill(Boundary,N)))
    @assert !isnothing(bnd_ind) throw(ArgumentError("must specify boundary"))
    bnd = find_type(Boundary,args...)
    simargs = (simargs..., bnd)

    dis_ind = findfirst(map(<:,types,fill(Discretization,N)))
    @assert !isnothing(dis_ind) throw(ArgumentError("must specify discretization"))
    dis = find_type(Discretization,args...)
    simargs = (simargs..., dis)

    sys = find_type(System,args...)
    sys = isnothing(sys) ? Scattering() : sys
    simargs = (simargs..., sys)

    sct = find_type(Scattering,args...)
    sct = isnothing(sct) ? Scattering() : sct
    simargs = (simargs..., sct)

    tls = find_type(TwoLevelSystem,args...)
    tls = isnothing(tls) ? TwoLevelSystem() : tls
    simargs = (simargs..., tls)

    lat = find_type(BravaisLattice,args...)
    lat = isnothing(lat) ? BravaisLattice(bnd) : lat
    simargs = (simargs..., lat)

    return Simulation(simargs...;kwargs...)
end

Simulation(sim::Simulation; sys=sim.sys, bnd=sim.bnd, dis=sim.dis, sct=sim.sct, tls=sim.tls, lat=sim.lat, disp_opt=false) = deepcopy(Simulation(sys=System(sys), bnd=Boundary(bnd), dis=Discretization(dis), sct=Scattering(sct), tls=TwoLevelSystem(tls), lat=BravaisLattice(lat), disp_opt=disp_opt))
