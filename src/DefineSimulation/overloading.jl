"""
    sys = System(domain1, domain2, ...)
"""
System(args::Vararg{Domain}) = System(args)


"""
    dis = Discretization{coordinate_system=Cartesian}(dx[, sub_pixel_num=DEFAULT_SUBSAMPLE_NUMBER])

discretization object for Simulation

`dx` lattice spacing (scalar)

`sub_pixel_num` is the subsampling rate used in sub-pixel smoothing.
"""
Discretization(dx::Union{Number,Array,Tuple}) = Discretization(dx, Cartesian())
Discretization(dx::Real, coordinate_system::CS, sub_pixel_num::Int=DEFAULT_SUBSAMPLE_NUMBER) where CS<:CoordinateSystem = Discretization((dx,dx), coordinate_system, sub_pixel_num, (0.,0.), (1,1), ((0,0),(0,0)))
Discretization(dx::Union{Array,Tuple}, coordinate_system::CS, sub_pixel_num::Int=DEFAULT_SUBSAMPLE_NUMBER) where CS<:CoordinateSystem = Discretization(Tuple(dx), coordinate_system, sub_pixel_num, (0.,0.), (1,1), ((0,0),(0,0)))

################################################################################
### BOUNDARY
################################################################################
Boundary(∂Ω, bcls...) = Boundary(_bnd_Ω(∂Ω),_bnd_bcls(bcls)...)

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

_find_bcs(bcls::Tuple) = _find_bcs(bcls...)
_find_bcs(bcl::AbstractBC,bcls...) = (bcl,_find_bcs(bcls...)...)
_find_bcs(bcl,bcls...) = (_find_bcs(bcls...)...,)
_find_bcs(bcl::AbstractBC) = (bcl,)
_find_bcs(bcl) = ()

_find_bls(bcls::Tuple) = _find_bls(bcls...)
_find_bls(bcl::AbstractBL,bcls...) = (bcl,_find_bls(bcls...)...)
_find_bls(bcl,bcls...) = (_find_bls(bcls...)...,)
_find_bls(bcl::AbstractBL) = (bcl,)
_find_bls(bcl) = ()

function _bnd_bcls(bcls::Vararg{Union{AbstractBC,AbstractBL},N}) where {N}
     bcs = _find_bcs(bcls)
     @assert length(bcs)==4 "insufficient number of boundaries ($(length(bcs))) specified"
     bcs = reorder(bcs)

     bls = _find_bls(bcls)
     @assert length(bls)==4 "insufficient number of boundaries layers ($(length(bls))) specified"
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
function Boundary(;∂Ω=fill(NaN,4), bc::TBC=DirichletBC(), bl::TBL=noBL(), depth::Real=0, r=false, kwargs...) where {TBC<:AbstractBC,TBL<:AbstractBL}
    if !isnan(r)
        @assert TBC<:AbstractBC "type of bc is $bc. when radius is specified, bc must be a boundary condition"
        @assert TBL<:AbstractBL "type of bl is $bl. when radius is specified, bl must be a boundary layer"
        ∂Ω = (0, r, 0, 2π)
        if TBC<:MatchedBC
            BC = MatchedBC{1,2}(kwargs[:outgoing_qns],kwargs[:incoming_qns])
        elseif TBC<:DirichletBC
            BC = DirichletBC{1,2}()
        elseif TBC<:NeumannBC
            BC = NeumannBC{1,2}()
        elseif TBC<:RobinBC
            BC = RobinBC{1,2}()
        else
            throw(Exception("polar boundaries not defined for BC $(TBC)"))
        end
        BCs = ( DirichletBC{1,1}(), BC, FloquetBC{2,1}(BravaisLattice(b=2π)), FloquetBC{2,2}(BravaisLattice(b=2π)) )
        if TBL<:noBL
            BL = noBL{1,2}(depth)
        elseif TBL<:PML
            BL = PML{1,2}(depth)
        elseif TBL<:cPML
            BL = cPML{1,2}(depth)
        else
            throw(Exception("unrecognized boundary layer $(TBL)"))
        end
        BLs = ( noBL{1,1}(0), BL, noBL{2,1}(0), noBL{2,2}(0) )
    else
        if TBC == DirichletBC
            BCs = (DirichletBC{1,1}(),DirichletBC{1,2}(),DirichletBC{2,1}(),DirichletBC{2,2}())
        elseif bc == NeumannBC
            BCs = (NeumannBC{1,1}(),NeumannBC{1,2}(),NeumannBC{2,1}(),NeumannBC{2,2}())
        elseif bc == FloquetBC
            BCs = (FloquetBC{1,1}(),FloquetBC{1,2}(),FloquetBC{2,1}(),FloquetBC{2,2}())
        else
            throw(ArgumentError("$(bc) incompatible with this wrapper signature of Boundary. For more control, use another signature"))
        end

        if TBL == PML
            BLs = (PML{1,1}(depth),PML{1,2}(depth),PML{2,1}(depth),PML{2,2}(depth))
        elseif TBL == cPML
            BLs = (cPML{1,1}(depth),cPML{1,2}(depth),cPML{2,1}(depth),cPML{2,2}(depth))
        elseif TBL == noBL
            BLs = (noBL{1,1}(),noBL{1,2}(),noBL{2,1}(),noBL{2,2}())
        else
            throw(ArgumentError("unrecognized bl $(bl)"))
        end
    end
    return Boundary(∂Ω, BCs..., BLs...)
end


################################################################################
### SCATTERING
################################################################################
Scattering(channel::Channels) = Scattering([channel])
"""
    sct = Scattering(channel1, channel2, ...)
"""
Scattering(args::Vararg{Channels}) = Scattering(vcat(args...))

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
    sys = isnothing(sys) ? System() : sys
    simargs = (simargs..., sys)

    sct = find_type(Scattering,args...)
    sct = isnothing(sct) ? Scattering() : sct
    simargs = (simargs..., sct)

    lat = find_type(BravaisLattice,args...)
    lat = isnothing(lat) ? BravaisLattice(bnd) : lat
    simargs = (simargs..., lat)

    return Simulation(simargs...;kwargs...)
end
