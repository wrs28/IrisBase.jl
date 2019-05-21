module BoundaryConditions

const EXTINCTION = 1e-10 # extinction in PML layer
const SCALING_ANGLE = .25 # phase in conductivity to accelerate evanscent decay

using ..Bravais,
..CoordinateSystems,
Formatting,
LinearAlgebra,
SparseArrays,
SpecialFunctions

export AbstractBL,
PML,
cPML,
noBL,
AbstractBC,
AbstractLocalBC,
DirichletBC,
NeumannBC,
RobinBC,
FloquetBC,
MatchedBC,
isDirichlet,
isNeumann,
isMatched,
isFloquet,
isPML,
iscPML,
isnoBL
# ezbc,
# ezbl


################################################################################
### BOUNDARY CONDITIONS
################################################################################
abstract type AbstractBL{DIM,SIDE} end
abstract type AbstractBC{DIM,SIDE} end
abstract type AbstractLocalBC{DIM,SIDE}<:AbstractBC{DIM,SIDE} end

################## BOUNDARY LAYERS #####################
"""
    PML{DIM,SIDE}(depth[,Σ])
"""
mutable struct PML{DIM,SIDE,FIN}<:AbstractBL{DIM,SIDE}
    depth::Float64 # depth of PML layer
    Δ::NTuple{2,NTuple{2,Float64}} # system size w/o pml
    bounds::Array{Float64,1}

    # create PML{DIM,SIDE} empty
    PML() = new{Union{},Union{},false}()
    PML{DIM,SIDE}() where {DIM,SIDE} = new{DIM,SIDE,false}(NaN,((NaN,NaN),(NaN,NaN)),[-Inf,Inf])
    # create PML{DIM,SIDE} initialized with depth
    PML{DIM,SIDE}(depth::Number,bounds::Array=[-Inf,Inf]) where {DIM,SIDE} = new{DIM,SIDE,true}(depth,((NaN,NaN),(NaN,NaN)),bounds)
    function PML{DIM,SIDE}(depth::Number,Δ::NTuple{2,NTuple{2,Number}},bounds::Array) where {DIM,SIDE}
        @assert depth>0 "cannot have PML with 0 depth"
        new{DIM,SIDE,true}(depth,Δ,bounds)
    end
    # add depth and and system size to PML
    (pml::PML{DIM,SIDE,false})(depth::Number,Δ::NTuple{2,NTuple{2,Number}},bounds::Array) where {DIM,SIDE} = PML{DIM,SIDE}(depth,Δ,bounds)
    # add system size to initialized PML even if given depth
    (pml::PML{DIM,SIDE,true})(depth::Number,Δ::NTuple{2,NTuple{2,Number}},bounds::Array) where {DIM,SIDE} = PML{DIM,SIDE}(pml.depth,Δ,pml.bounds)
    # add system size to initialized PML
    (pml::PML{DIM,SIDE,true})(Δ::NTuple{2,NTuple{2,Number}}) where {DIM,SIDE} = PML{DIM,SIDE}(pml.depth,Δ,pml.bounds)
    # if PML initialized, do nothing unless given depth and/or system size
    (pml::PML{DIM,SIDE,true})(args...) where {DIM,SIDE} = pml

    (pml::PML{DIM,SIDE,true})(x::Number) where {DIM,SIDE} = conductivity_profile(x,pml.Δ,pml.depth,DIM,SIDE)
    function (pml::PML{1,SIDE,true})(x::Number,y::Number) where {SIDE}
        if pml.bounds[1] ≤ y ≤ pml.bounds[2]
            return conductivity_profile(x,pml.Δ,pml.depth,1,SIDE)
        else
            return conj_conductivity_profile(x,pml.Δ,pml.depth,1,SIDE)
        end
    end
    function (pml::PML{2,SIDE,true})(x::Number,y::Number) where {SIDE}
        if pml.bounds[1] ≤ x ≤ pml.bounds[2]
            return conductivity_profile(y,pml.Δ,pml.depth,2,SIDE)
        else
            return conj_conductivity_profile(y,pml.Δ,pml.depth,2,SIDE)
        end
    end

    Base.show(io::IO,::PML{Union{}}) = print("PML()")
    Base.show(io::IO,pml::PML{DIM,SIDE}) where {DIM,SIDE} = print("PML{$DIM,$SIDE}($(pml.depth))")
end
"""
    cPML{DIM,SIDE}(depth[,Σ])
"""
mutable struct cPML{DIM,SIDE,FIN}<:AbstractBL{DIM,SIDE}
    depth::Float64
    Δ::NTuple{2,NTuple{2,Float64}} # tuple of tuple of depths
    bounds::Array{Float64,1}

    cPML() = new{Union{},Union{},false}()
    cPML{DIM,SIDE}() where {DIM,SIDE} = new{DIM,SIDE,false}(NaN,((NaN,NaN),(NaN,NaN)),[-Inf,Inf])
    cPML{DIM,SIDE}(depth::Number,bounds::Array=[-Inf,Inf]) where {DIM,SIDE} = new{DIM,SIDE,true}(depth,((NaN,NaN),(NaN,NaN)),bounds)
    function cPML{DIM,SIDE}(depth::Number,Δ::NTuple{2,NTuple{2,Number}},bounds::Array) where {DIM,SIDE}
        @assert depth>0 "cannot have cPML with 0 depth"
        new{DIM,SIDE,true}(depth,Δ,bounds)
    end
    # add depth and and system size to PML
    (cpml::cPML{DIM,SIDE,false})(depth::Number,Δ::NTuple{2,NTuple{2,Number}},bounds::Array) where {DIM,SIDE} = cPML{DIM,SIDE}(depth,Δ,bounds)
    # add system size to initialized PML even if given depth
    (cpml::cPML{DIM,SIDE,true})(depth::Number,Δ::NTuple{2,NTuple{2,Number}},bounds::Array) where {DIM,SIDE} = cPML{DIM,SIDE}(cpml.depth,Δ,cpml.bounds)
    # add system size to initialized PML
    (cpml::cPML{DIM,SIDE,true})(Δ::NTuple{2,NTuple{2,Number}}) where {DIM,SIDE} = cPML{DIM,SIDE}(cpml.depth,Δ,cpml.bounds)
    # if PML initialized, do nothing unless given depth and/or system size
    (cpml::cPML{DIM,SIDE,true})(args...) where {DIM,SIDE} = cpml

    (cpml::cPML{DIM,SIDE,true})(x::Number) where {DIM,SIDE} = conj_conductivity_profile(x,cpml.Δ,cpml.depth,DIM,SIDE)
    function (cpml::cPML{1,SIDE,true})(x::Number,y::Number) where {SIDE}
        if cpml.bounds[1] ≤ y ≤ cpml.bounds[2]
            return conj_conductivity_profile(x,cpml.Δ,cpml.depth,1,SIDE)
        else
            return conductivity_profile(x,cpml.Δ,cpml.depth,1,SIDE)
        end
    end
    function (cpml::cPML{2,SIDE,true})(x::Number,y::Number) where {SIDE}
        if cpml.bounds[1] ≤ x ≤ cpml.bounds[2]
            return conj_conductivity_profile(y,cpml.Δ,cpml.depth,2,SIDE)
        else
            return conductivity_profile(y,cpml.Δ,cpml.depth,2,SIDE)
        end
    end

    Base.show(io::IO,::cPML{Union{}}) = print("cPML()")
    Base.show(io::IO,cpml::cPML{DIM,SIDE}) where {DIM,SIDE} = print("cPML{$DIM,$SIDE}($(cpml.depth))")
end
"""
    noBL{DIM,SIDE}()
"""
mutable struct noBL{DIM,SIDE}<:AbstractBL{DIM,SIDE}
    depth::Float64

    noBL() = new{Union{},Union{}}()
    noBL{DIM,SIDE}(args...) where {DIM,SIDE} = new{DIM,SIDE}(0)
    (nbl::noBL)(depth::Number,Δ::Tuple) = nbl
    (nbl::noBL)(Δ::Tuple) = nbl
    (::noBL)(x::Number,args...) = 0
    Base.show(io::IO,::noBL{Union{}}) = print("noBL()")
    Base.show(io::IO,cpml::noBL{DIM,SIDE}) where {DIM,SIDE} = print("noBL{$DIM,$SIDE}")
end


function conductivity_profile(x,Δ,depth,dim,side)
    β = 2*sqrt(depth)
    α = -(1/4)*exp(complex(0,SCALING_ANGLE))*log(EXTINCTION)/depth/((2exp(β)-β)/(2β*expm1(β))-1/β^2)
    D = Δ[dim][side]-depth*(-1)^side
    if sign(D-x)*(-1)^side ≤ 0
        u = abs(x-D)/depth
        σ = α*u*expm1(β*u)/expm1(β)
    else
        σ = complex(0.)
    end
    #BL[i][j]<:cPML ? α[k] = flipsign(conj(α[k]),-1) : nothing
    return isnan(σ) ? complex(0.) : σ
end
conj_conductivity_profile(x,Δ,depth,DIM,SIDE) = -conj(conductivity_profile(x,Δ,depth,DIM,SIDE))

################## BOUNDARY CONDITIONS #####################
"""
    dbc = DirichletBC{DIM,SIDE}(N[,g])
"""
struct DirichletBC{DIM,SIDE,FIN,T,Tg,Tgf,TS} <: AbstractLocalBC{DIM,SIDE}
    g::Tg
    gfun::Tgf
    N::NTuple{2,Int}
    xmin::NTuple{2,Float64}
    I1::Array{Array{Int,1},1}
    J1::Array{Array{Int,1},1}
    I2::Array{Array{Int,1},1}
    J2::Array{Array{Int,1},1}
    weights::Array{Array{Float64,1},1}
    S::TS

    DirichletBC() = new{Union{},Union{},false,Union{},Union{},Union{},Union{}}()
    DirichletBC{DIM,SIDE}() where {DIM,SIDE} = new{DIM,SIDE,false,Float64,Float64,typeof(identity),Union{}}(0.0,identity,(0,0),(0.0,0.0))
    DirichletBC{DIM,SIDE}(g::Tg) where {DIM,SIDE,Tg} = new{DIM,SIDE,false,Tg,Tg,typeof(identity),Union{}}(g,identity,(0,0),(0.0,0.0))
    DirichletBC{DIM,SIDE}(g::Function) where {DIM,SIDE} = new{DIM,SIDE,false,typeof(g),Float64,typeof(g),Union{}}(NaN,g,(0,0),(0.0,0.0))
    function DirichletBC{DIM,SIDE}(g::Number,N::Union{Array,Tuple}; kwargs...) where {DIM,SIDE}
        dim_perp = mod1(DIM+1,2)
        N_perp = N[dim_perp]
        typeof(g)<:Number ? g = fill(g,N_perp) : nothing
        return DirichletBC{DIM,SIDE}(g,N; kwargs...)
    end
    function DirichletBC{DIM,SIDE}(g::Function,N::Union{Array,Tuple}; kwargs...) where {DIM,SIDE}
        xmin = haskey(kwargs,:xmin) ? kwargs[:xmin] : (0.0,0.0)
        @assert haskey(kwargs,:dx) "DirichletBC(function,...) requires kwarg :dx"
        dx = kwargs[:dx]
        dim_perp = mod1(DIM+1,2)
        x = xmin[dim_perp] .+ ((1:N[dim_perp]) .- 1/2)*dx[dim_perp]
        DirichletBC{DIM,SIDE}(g.(x),N; kwargs...)
    end
    function DirichletBC{DIM,SIDE}(g::Array,N::Union{Array,Tuple}; kwargs...) where {DIM,SIDE}
        @assert DIM ≤ 2 "unrecognized dimension d=$DIM. must have d≤2 for now."
        xmin = haskey(kwargs,:xmin) ? kwargs[:xmin] : (0.0,0.0)
        N_perp = N[mod1(DIM+1,2)]
        I1, I2 = [Array{Int}(undef,N_perp)], [Array{Int}(undef,N_perp)]
        J1, J2 = [Array{Int}(undef,N_perp)], [Array{Int}(undef,N_perp)]
        weights = [Array{eltype(g)}(undef,N_perp)]
        for idx ∈ 1:N_perp
            if DIM==1
                I1[1][idx] = SIDE==1 ? 1 : N[1]
                J1[1][idx] = idx
                I2[1][idx] = SIDE==1 ? 1 : N[1]
                J2[1][idx] = idx
            else
                I1[1][idx] = idx
                J1[1][idx] = SIDE==1 ? 1 : N[2]
                I2[1][idx] = idx
                J2[1][idx] = SIDE==1 ? 1 : N[2]
            end
            weights[1][idx] = -1
        end
        S = 2g
        return new{DIM,SIDE,true,typeof(g),typeof(g),typeof(identity),typeof(S)}(g,identity,Tuple(N),Tuple(xmin),I1,J1,I2,J2,weights,S)
    end

    (dbc::DirichletBC{DIM,SIDE,false})(N::Union{Tuple,Array}; kwargs...) where {DIM,SIDE} = DirichletBC{DIM,SIDE}(dbc.g, N; kwargs...)
    (dbc::DirichletBC{DIM,SIDE,Tg,false})(N::Union{Tuple,Array}; kwargs...) where {DIM,SIDE,Tg<:Function} = DirichletBC{DIM,SIDE}(dbc.gfun, N; kwargs...)
    (dbc::DirichletBC{DIM,SIDE,true})(args...; kwargs...) where {DIM,SIDE} = dbc

    (dbc::DirichletBC{DIM,SIDE,true})(k::Number,args...) where {DIM,SIDE} = [1]

    Base.show(io::IO,::DirichletBC{Union{},Union{},false}) = print("DirichletBC()")
    Base.show(io::IO,DBC::DirichletBC{DIM,SIDE,true}) where {DIM,SIDE} = print("DirichletBC{$DIM,$SIDE}")
    Base.show(io::IO,DBC::DirichletBC{DIM,SIDE,false}) where {DIM,SIDE} = print("DirichletBC{$DIM,$SIDE}")
    Base.show(io::IO,DBC::DirichletBC{DIM,SIDE,false,T}) where {DIM,SIDE,T<:Number} = print("DirichletBC{$DIM,$SIDE}(",fmt("2.1f",DBC.g),")")
    Base.show(io::IO,DBC::DirichletBC{DIM,SIDE,false,T}) where {DIM,SIDE,T<:Function} = print("DirichletBC{$DIM,$SIDE}($(DBC.gfun))")
end

"""
    nbc = NeumannBC{DIM,SIDE}(N[,g])
"""
struct NeumannBC{DIM,SIDE,FIN,T,Tg,Tgf,TS} <: AbstractLocalBC{DIM,SIDE}
    g::Tg
    gfun::Tgf
    N::NTuple{2,Int}
    xmin::NTuple{2,Float64}
    I1::Array{Array{Int,1},1}
    J1::Array{Array{Int,1},1}
    I2::Array{Array{Int,1},1}
    J2::Array{Array{Int,1},1}
    weights::Array{Array{Float64,1},1}
    S::TS

    NeumannBC() = new{Union{},Union{},false,Union{},Union{},Union{},Union{}}()
    NeumannBC{DIM,SIDE}() where {DIM,SIDE} = new{DIM,SIDE,false,Float64,Float64,typeof(identity),Union{}}(0.0,identity,(0,0),(0.0,0.0))
    NeumannBC{DIM,SIDE}(g::Tg) where {DIM,SIDE,Tg} = new{DIM,SIDE,false,Tg,Tg,typeof(identity),Union{}}(g,identity,(0,0),(0.0,0.0))
    NeumannBC{DIM,SIDE}(g::Function) where {DIM,SIDE} = new{DIM,SIDE,false,typeof(g),Float64,typeof(g),Union{}}(NaN,g,(0,0),(0.0,0.0))
    function NeumannBC{DIM,SIDE}(g::Number,N::Union{Array,Tuple}; kwargs...) where {DIM,SIDE}
        dim_perp = mod1(DIM+1,2)
        N_perp = N[dim_perp]
        typeof(g)<:Number ? g = fill(g,N_perp) : nothing
        return NeumannBC{DIM,SIDE}(g,N; kwargs...)
    end
    function NeumannBC{DIM,SIDE}(g::Function,N::Union{Array,Tuple}; kwargs...) where {DIM,SIDE}
        xmin = haskey(kwargs,:xmin) ? kwargs[:xmin] : (0.0,0.0)
        @assert haskey(kwargs,:dx) "NeumannBC requires kwarg :dx"
        dx = kwargs[:dx]
        dim_perp = mod1(DIM+1,2)
        x = xmin[dim_perp] .+ ((1:N[dim_perp]) .- 1/2)*dx[dim_perp]
        return NeumannBC{DIM,SIDE}(g.(x),N; kwargs...)
    end
    function NeumannBC{DIM,SIDE}(g::Array,N::Union{Array,Tuple}; kwargs...) where {DIM,SIDE}
        @assert DIM ≤ 2 "unrecognized dimension d=$DIM. must have d≤2 for now."
        @assert haskey(kwargs,:dx) "NeumannBC requires kwarg :dx"
        dx = kwargs[:dx]
        xmin = haskey(kwargs,:xmin) ? kwargs[:xmin] : (0.0,0.0)
        N_perp = N[mod1(DIM+1,2)]
        I1, I2 = [Array{Int}(undef,N_perp)], [Array{Int}(undef,N_perp)]
        J1, J2 = [Array{Int}(undef,N_perp)], [Array{Int}(undef,N_perp)]
        weights = [Array{eltype(g)}(undef,N_perp)]
        for idx ∈ 1:N_perp
            if DIM==1
                I1[1][idx] = SIDE==1 ? 1 : N[1]
                J1[1][idx] = idx
                I2[1][idx] = SIDE==1 ? 1 : N[1]
                J2[1][idx] = idx
            else
                I1[1][idx] = idx
                J1[1][idx] = SIDE==1 ? 1 : N[2]
                I2[1][idx] = idx
                J2[1][idx] = SIDE==1 ? 1 : N[2]
            end
            weights[1][idx] = 1
        end
        S = 2g/dx[DIM]
        return new{DIM,SIDE,true,typeof(g),typeof(g),typeof(identity),typeof(S)}(g,identity,Tuple(N),Tuple(xmin),I1,J1,I2,J2,weights,S)
    end

    (nbc::NeumannBC{DIM,SIDE,false})(N::Union{Tuple,Array}; kwargs...) where {DIM,SIDE} = NeumannBC{DIM,SIDE}(nbc.g, N; kwargs...)
    (nbc::NeumannBC{DIM,SIDE,false,Tg})(N::Union{Tuple,Array}; kwargs...) where {DIM,SIDE,Tg<:Function} = NeumannBC{DIM,SIDE}(nbc.gfun, N; kwargs...)
    (nbc::NeumannBC{DIM,SIDE,true})(args...; kwargs...) where {DIM,SIDE} = nbc

    (nbc::NeumannBC{DIM,SIDE,true})(k::Number,args...) where {DIM,SIDE} = [1]

    Base.show(io::IO,::NeumannBC{Union{},Union{},false}) = print("NeumannBC()")
    Base.show(io::IO,::NeumannBC{DIM,SIDE,true}) where {DIM,SIDE} = print("NeumannBC{$DIM,$SIDE}")
    Base.show(io::IO,::NeumannBC{DIM,SIDE,false}) where {DIM,SIDE} = print("NeumannBC{$DIM,$SIDE}")
    Base.show(io::IO,NBC::NeumannBC{DIM,SIDE,false,T}) where {DIM,SIDE,T<:Number} = print("NeumannBC{$DIM,$SIDE}(",fmt("2.1f",NBC.g),")")
    Base.show(io::IO,NBC::NeumannBC{DIM,SIDE,false,T}) where {DIM,SIDE,T<:Function} = print("NeumannBC{$DIM,$SIDE}($(NBC.gfun))")
end

"""
    rbc = RobinBC{DIM,SIDE}(N, α, β, dx[, g])

Robin Boundary conditions A*ϕ+B*∇ϕ = G, all functions of position
"""
struct RobinBC{DIM,SIDE,FIN,TA,TB,TG,Ta,Tb,Tg,Taf,Tbf,Tgf,TW,TS} <: AbstractLocalBC{DIM,SIDE}
    α::Ta
    β::Tb
    g::Tg

    αfun::Taf
    βfun::Tbf
    gfun::Tgf

    N::NTuple{2,Int}
    xmin::NTuple{2,Float64}
    I1::Array{Array{Int,1},1}
    J1::Array{Array{Int,1},1}
    I2::Array{Array{Int,1},1}
    J2::Array{Array{Int,1},1}
    weights::Array{Array{TW,1},1}
    S::TS

    RobinBC() = new{Union{},Union{},false,Union{},Union{},Union{},Union{},Union{},Union{},Union{},Union{},Union{},Union{},Union{}}()
    RobinBC{DIM,SIDE}(α::Number,β::Number) where {DIM,SIDE} = RobinBC{DIM,SIDE}(α,β,zero(α*β))
    RobinBC{DIM,SIDE}(α::Function,β::Number) where {DIM,SIDE} = RobinBC{DIM,SIDE}(α,β,zero(α(1.)*β))
    RobinBC{DIM,SIDE}(α::Number,β::Function) where {DIM,SIDE} = RobinBC{DIM,SIDE}(α,β,zero(α*β(1.)))
    RobinBC{DIM,SIDE}(α::Function,β::Function) where {DIM,SIDE} = RobinBC{DIM,SIDE}(α,β,zero(α(1.)*β(1.)))
    function RobinBC{DIM,SIDE}(α::Ta,β::Tb,g::Tg) where {DIM,SIDE,Ta,Tb,Tg}
        Ta<:Function ? (aval,afun) = (NaN,α) : (aval,afun) = (α,identity)
        Tb<:Function ? (bval,bfun) = (NaN,β) : (bval,bfun) = (β,identity)
        Tg<:Function ? (gval,gfun) = (NaN,g) : (gval,gfun) = (g,identity)
        new{DIM,SIDE,false,Ta,Tb,Tg,typeof(aval),typeof(bval),typeof(gval),typeof(afun),typeof(bfun),typeof(gfun),Union{},Union{}}(aval,bval,gval,afun,bfun,gfun,(0,0),(0.0,0.0))
    end
    function RobinBC{DIM,SIDE}(α,β,g,N::Union{Tuple,Array}; kwargs...) where {DIM,SIDE}
        return RobinBC{DIM,SIDE}(_αβg(α,DIM,N;kwargs...),_αβg(β,DIM,N;kwargs...),_αβg(g,DIM,N;kwargs...), N; kwargs...)
    end
    function RobinBC{DIM,SIDE}(α::AbstractArray, β::AbstractArray, g::AbstractArray, N::Union{Tuple,Array}; kwargs...) where {DIM,SIDE}
        @assert DIM ≤ 2 "unrecognized dimension d=$DIM. must have d≤2 for now."
        @assert haskey(kwargs,:dx) "RobinBC requires kwarg :dx"
        dx = kwargs[:dx]
        xmin = haskey(kwargs,:xmin) ? kwargs[:xmin] : (0.0,0.0)

        α, β, g = promote(α, β, g)
        A, B, G = Diagonal(α), Diagonal(β), g

        M_den = factorize(B+A*dx[DIM])
        M = M_den\(B-A*dx[DIM])
        S = +M_den\2G

        I1 = [Array{Int}(undef,prod(size(M)))]
        I2 = [Array{Int}(undef,prod(size(M)))]
        J1 = [Array{Int}(undef,prod(size(M)))]
        J2 = [Array{Int}(undef,prod(size(M)))]
        weights = [Array{eltype(M)}(undef,prod(size(M)))]
        idx = 1
        for k ∈ CartesianIndices(M)
            i,j = Tuple(k)
            if DIM==1
                I1[1][idx] = SIDE==1 ? 1 : N[1]
                J1[1][idx] = j
                I2[1][idx] = SIDE==1 ? 1 : N[1]
                J2[1][idx] = j
            else
                I1[1][idx] = i
                J1[1][idx] = SIDE==1 ? 1 : N[2]
                I2[1][idx] = i
                J2[1][idx] = SIDE==1 ? 1 : N[2]
            end
            weights[1][idx] = M[k]
            idx += 1
        end
        return new{DIM,SIDE,true,typeof(α),typeof(β),typeof(g),typeof(α),typeof(β),typeof(g),typeof(identity),typeof(identity),typeof(identity),eltype(weights[1]),typeof(S)}(α,β,g,identity,identity,identity,Tuple(N),Tuple(xmin),I1,J1,I2,J2,weights,S)
    end

    (rbc::RobinBC{DIM,SIDE,false,TA,TB,TG})(N; kwargs...) where {DIM,SIDE,TB,TG} where TA<:Function = RobinBC{DIM,SIDE}(rbc.αfun,rbc.β,rbc.g,N; kwargs...)
    (rbc::RobinBC{DIM,SIDE,false,TA,TB,TG})(N; kwargs...) where {DIM,SIDE,TA,TG} where TB<:Function = RobinBC{DIM,SIDE}(rbc.α,rbc.βfun,rbc.g,N; kwargs...)
    (rbc::RobinBC{DIM,SIDE,false,TA,TB,TG})(N; kwargs...) where {DIM,SIDE,TA,TB} where TG<:Function = RobinBC{DIM,SIDE}(rbc.α,rbc.β,rbc.gfun,N; kwargs...)
    (rbc::RobinBC{DIM,SIDE,false,TA,TB,TG})(N; kwargs...) where {DIM,SIDE,TG} where {TA<:Function,TB<:Function} = RobinBC{DIM,SIDE}(rbc.αfun,rbc.βfun,rbc.g,N; kwargs...)
    (rbc::RobinBC{DIM,SIDE,false,TA,TB,TG})(N; kwargs...) where {DIM,SIDE,TB} where {TA<:Function,TG<:Function} = RobinBC{DIM,SIDE}(rbc.αfun,rbc.β,rbc.gfun,N; kwargs...)
    (rbc::RobinBC{DIM,SIDE,false,TA,TB,TG})(N; kwargs...) where {DIM,SIDE,TA} where {TB<:Function,TG<:Function} = RobinBC{DIM,SIDE}(rbc.α,rbc.βfun,rbc.gfun,N; kwargs...)
    (rbc::RobinBC{DIM,SIDE,false,TA,TB,TG})(N; kwargs...) where {DIM,SIDE} where {TA<:Function,TB<:Function,TG<:Function} = RobinBC{DIM,SIDE}(rbc.αfun,rbc.βfun,rbc.gfun,N; kwargs...)
    (rbc::RobinBC{DIM,SIDE,false})(N; kwargs...) where {DIM,SIDE} = RobinBC{DIM,SIDE}(rbc.α,rbc.β,rbc.g,N; kwargs...)
    (rbc::RobinBC{DIM,SIDE,true})(args...; kwargs...) where {DIM,SIDE} = rbc

    (rbc::RobinBC{DIM,SIDE,true})(k::Number,args...) where {DIM,SIDE} = [1]

    Base.show(io::IO,::RobinBC{Union{},Union{},false,Union{},Union{},Union{},Union{}}) = print("RobinBC()")
    Base.show(io::IO,::RobinBC{DIM,SIDE,true}) where {DIM,SIDE} = print("RobinBC{$DIM,$SIDE}")
    Base.show(io::IO,::RobinBC{DIM,SIDE,false}) where {DIM,SIDE} = print("RobinBC{$DIM,$SIDE}")
    Base.show(io::IO,RBC::RobinBC{DIM,SIDE,false,TA,TB,TG}) where {DIM,SIDE} where {TA<:Number,TB<:Number,TG<:Number} = print("RobinBC{$DIM,$SIDE}($(RBC.α),$(RBC.β),$(RBC.g))")
    Base.show(io::IO,RBC::RobinBC{DIM,SIDE,false,TA,TB,TG}) where {DIM,SIDE} where {TA<:Function,TB<:Number,TG<:Number} = print("RobinBC{$DIM,$SIDE}($(RBC.αfun),$(RBC.β),$(RBC.g))")
    Base.show(io::IO,RBC::RobinBC{DIM,SIDE,false,TA,TB,TG}) where {DIM,SIDE} where {TA<:Number,TB<:Function,TG<:Number} = print("RobinBC{$DIM,$SIDE}($(RBC.α),$(RBC.βfun),$(RBC.g))")
    Base.show(io::IO,RBC::RobinBC{DIM,SIDE,false,TA,TB,TG}) where {DIM,SIDE} where {TA<:Number,TB<:Number,TG<:Function} = print("RobinBC{$DIM,$SIDE}($(RBC.α),$(RBC.β),$(RBC.gfun))")
    Base.show(io::IO,RBC::RobinBC{DIM,SIDE,false,TA,TB,TG}) where {DIM,SIDE} where {TA<:Function,TB<:Function,TG<:Number} = print("RobinBC{$DIM,$SIDE}($(RBC.αfun),$(RBC.βfun),$(RBC.g))")
    Base.show(io::IO,RBC::RobinBC{DIM,SIDE,false,TA,TB,TG}) where {DIM,SIDE} where {TA<:Function,TB<:Number,TG<:Function} = print("RobinBC{$DIM,$SIDE}($(RBC.αfun),$(RBC.β),$(RBC.gfun))")
    Base.show(io::IO,RBC::RobinBC{DIM,SIDE,false,TA,TB,TG}) where {DIM,SIDE} where {TA<:Number,TB<:Function,TG<:Function} = print("RobinBC{$DIM,$SIDE}($(RBC.α),$(RBC.βfun),$(RBC.gfun))")
    Base.show(io::IO,RBC::RobinBC{DIM,SIDE,false,TA,TB,TG}) where {DIM,SIDE} where {TA<:Function,TB<:Function,TG<:Function} = print("RobinBC{$DIM,$SIDE}($(RBC.αfun),$(RBC.βfun),$(RBC.gfun))")
end
function _αβg(f::T,DIM,N;kwargs...) where T<:Function
    @assert haskey(kwargs,:dx) "need to give dx as keyword"
    xmin = haskey(kwargs,:xmin) ? kwargs[:xmin] : (0.0,0.0)
    @assert haskey(kwargs,:dx) "RobinBC requires kwarg :dx"
    dx = kwargs[:dx]
    dim_perp = mod1(DIM+1,2)
    x = xmin[dim_perp] .+ ((1:N[dim_perp]) .- 1/2)*dx[dim_perp]
    return f.(x)
end
function _αβg(f::Number,DIM,N;kwargs...)
    dim_perp = mod1(DIM+1,2)
    N_perp = N[dim_perp]
    return fill(f,N_perp)
end
_αβg(x) = x


"""
    pbc = FloquetBC{DIM,SIDE}(N[, lattice=BravaisLattice(a=N[1],b=N[2])])
"""
struct FloquetBC{DIM,SIDE,FIN} <: AbstractBC{DIM,SIDE}
    lattice::BravaisLattice
    N::NTuple{2,Int}
    I1::Array{Array{Array{Int,1},1},1}
    J1::Array{Array{Array{Int,1},1},1}
    I2::Array{Array{Array{Int,1},1},1}
    J2::Array{Array{Array{Int,1},1},1}
    weights::Array{Array{Array{Float64,1},1},1}
    shifts::Array{Array{Int,1},1}

    FloquetBC() = new{Union{},Union{},false}()
    FloquetBC{DIM,SIDE}() where {DIM,SIDE} = new{DIM,SIDE,false}()
    FloquetBC{DIM,SIDE}(lattice::BravaisLattice; kwargs...) where {DIM,SIDE} = new{DIM,SIDE,false}(lattice,(0,0))
    function FloquetBC{DIM,SIDE}(lattice,N; kwargs...) where {DIM,SIDE}

        @assert !(0 ∈ N) "improperly initialized N=$N"

        a, α, v1 = lattice.a, lattice.α, lattice.v1
        b, β, v2 = lattice.b, lattice.β, lattice.v2
        width, height  = a*sin(β-α)/sin(β), b*sin(β)
        dx, dy = width/N[1], height/N[2]

        @assert DIM ∈ [1,2] "invalid dimeions $(DIM)"
        if (DIM == 1 && isinf(a)) || (DIM == 2 && isinf(b))
            return new{DIM,SIDE,true}(lattice, N, [Int[]], [Int[]], [Int[]], [Int[]], [Float64[]], [Int[]])
        end

        if !isinf(a)
            start, stop = lattice.x0 + min(0, a*sin(β-α)/sin(β)), lattice.x0 + max(0, a*sin(β-α)/sin(β))
            x = LinRange(start + dx/2, stop - dx/2, N[1])
        else
            x = 1/2:N[1]-1/2#LinRange(lattice.x0,lattice.x0,1)
            dx = 1
        end

        if !isinf(b)
            start, stop = lattice.y0 + min(0, b*sin(β)), lattice.y0 + max(0, b*sin(β))
            y = LinRange(start + dy/2, stop - dy/2, N[2])
        else
            y = 1/2:N[2]-1/2#LinRange(lattice.y0,lattice.y0,1)
            dy = 1
        end

        @assert SIDE ∈ [1,2] "unrecognized side $side, must be 1 or 2"
        if DIM == 1
            p1, p2 = SIDE==1 ? bravais_coordinates(x[1]-dx, y, lattice) : bravais_coordinates(x[end]+dx, y, lattice)
        else
            p1, p2 = SIDE==1 ? bravais_coordinates(x, y[1]-dy, lattice) : bravais_coordinates(x, y[end]+dy, lattice)
        end
        Ma, Mb = -floor.(Int, p1/a), -floor.(Int, p2/b)

        if isinf(a)
            X = v1[1]*p1 + v2[1]*(p2 + Mb*b)
            Y = v1[2]*p1 + v2[2]*(p2 + Mb*b)
        elseif isinf(b)
            X = v1[1]*(p1 + Ma*a) + v2[1]*p2
            Y = v1[2]*(p1 + Ma*a) + v2[2]*p2
        else
            X = v1[1]*(p1 + Ma*a) + v2[1]*(p2 + Mb*b)
            Y = v1[2]*(p1 + Ma*a) + v2[2]*(p2 + Mb*b)
        end

        Ma += -floor.(Int, X/width)
        Mb += -floor.(Int, Y/height)

        x_inds1 = floor.(Int, X/dx .+ 1/2)
        x_inds2 = x_inds1 .+ 1
        y_inds1 = floor.(Int, Y/dy .+ 1/2)
        y_inds2 = y_inds1 .+ 1

        Cx1, Cx2 = abs.(X/dx .+ 1/2 - x_inds2), abs.(X/dx .+ 1/2 - x_inds1)
        cx1, cx2 = Cx1./(Cx1+Cx2), Cx2./(Cx1+Cx2)

        Cy1, Cy2 = abs.(Y/dy .+ 1/2 - y_inds2), abs.(Y/dy .+ 1/2 - y_inds1)
        cy1, cy2 = Cy1./(Cy1+Cy2), Cy2./(Cy1+Cy2)

        q, r, s = Array{Int}(undef,2), Array{Int}(undef,2), Array{Float64}(undef,2)
        t, u, v = Array{Int}(undef,2), Array{Int}(undef,2), Array{Float64}(undef,2)

        I1 = Array{Int}(undef, 4N[mod1(DIM+1,2)])
        I2 = Array{Int}(undef, 4N[mod1(DIM+1,2)])
        J1 = Array{Int}(undef, 4N[mod1(DIM+1,2)])
        J2 = Array{Int}(undef, 4N[mod1(DIM+1,2)])
        V = Array{Float64}(undef, 4N[mod1(DIM+1,2)])
        na = Array{Int}(undef, 4N[mod1(DIM+1,2)])
        nb = Array{Int}(undef, 4N[mod1(DIM+1,2)])

        for i ∈ 1:N[mod1(DIM+1,2)]

            if SIDE==1
                (ind_x, ind_y) = DIM==1 ? (1, i) : (i, 1)
            else
                (ind_x, ind_y) = DIM==1 ? (N[1], i) : (i, N[2])
            end

            q[1], q[2] = ind_x, ind_x
            r[1], r[2] = mod1(x_inds1[i],N[1]), mod1(x_inds2[i],N[1])
            s[1], s[2] = cx1[i], cx2[i]

            t[1], t[2] = ind_y, ind_y
            u[1], u[2] = mod1(y_inds1[i],N[2]), mod1(y_inds2[i],N[2])
            v[1], v[2] = cy1[i], cy2[i]

            idx = 1
            for j ∈ 1:2, k ∈ 1:2
                I1[4(i-1)+idx] = q[j]
                I2[4(i-1)+idx] = r[j]
                J1[4(i-1)+idx] = t[k]
                J2[4(i-1)+idx] = u[k]
                V[4(i-1)+idx] = s[j]*v[k]
                idx += 1
            end
            na[(4(i-1)+1):(4(i-1)+4)] .= Ma[i]
            nb[(4(i-1)+1):(4(i-1)+4)] .= Mb[i]
        end
        Na = findall.(isequal.(unique(na)),[na])
        Nb = findall.(isequal.(unique(nb)),[nb])
        Naa = Array{Int}(undef,length(Na))
        Nbb = Array{Int}(undef,length(Nb))
        for i ∈ 1:length(Na)-1
            Naa[i] = na[Na[i][1]]
        end
        for i ∈ 1:length(Nb)-1
            Nbb[i] = nb[Nb[i][1]]
        end
        Naa[end] = na[Na[end][1]]
        Nbb[end] = nb[Nb[end][1]]
        i1 = Array{Array{Array{Int,1},1},1}(undef,2)
        i2 = Array{Array{Array{Int,1},1},1}(undef,2)
        j1 = Array{Array{Array{Int,1},1},1}(undef,2)
        j2 = Array{Array{Array{Int,1},1},1}(undef,2)
        v = Array{Array{Array{Float64,1},1},1}(undef,2)
        i1[1] = Array{Array{Int,1},1}(undef,length(Na))
        i2[1] = Array{Array{Int,1},1}(undef,length(Na))
        j1[1] = Array{Array{Int,1},1}(undef,length(Na))
        j2[1] = Array{Array{Int,1},1}(undef,length(Na))
        v[1] = Array{Array{Float64,1},1}(undef,length(Na))
        for n ∈ eachindex(Na)
            i1[1][n] = I1[Na[n]]
            i2[1][n] = I2[Na[n]]
            j1[1][n] = J1[Na[n]]
            j2[1][n] = J2[Na[n]]
            v[1][n] = V[Na[n]]
        end
        i1[2] = Array{Array{Int,1},1}(undef,length(Nb))
        i2[2] = Array{Array{Int,1},1}(undef,length(Nb))
        j1[2] = Array{Array{Int,1},1}(undef,length(Nb))
        j2[2] = Array{Array{Int,1},1}(undef,length(Nb))
        v[2] = Array{Array{Float64,1},1}(undef,length(Nb))
        for n ∈ eachindex(Nb)
            i1[2][n] = I1[Nb[n]]
            i2[2][n] = I2[Nb[n]]
            j1[2][n] = J1[Nb[n]]
            j2[2][n] = J2[Nb[n]]
            v[2][n] = V[Nb[n]]
        end
        return new{DIM,SIDE,true}(lattice, Tuple(N), i1, j1, i2, j2, v, [Naa, Nbb])
    end

    (pbc::FloquetBC{DIM,SIDE,false})(N; kwargs...) where {DIM,SIDE} = FloquetBC{DIM,SIDE}(pbc.lattice,N; kwargs...)
    (pbc::FloquetBC{DIM,SIDE,false})(lattice::BravaisLattice, N; kwargs...) where {DIM,SIDE} = FloquetBC{DIM,SIDE}(lattice,N; kwargs...)
    (pbc::FloquetBC{DIM,SIDE,true})(args...;kwargs...) where {DIM,SIDE} = pbc

    function (pbc::FloquetBC{DIM,SIDE,true})(k::Number,ka::Number,kb::Number,ind::Int) where {DIM,SIDE}
        (K,k) = ind==1 ? (ka,kb) : (kb,ka)
        (R,r) = ind==1 ? (pbc.lattice.a,pbc.lattice.b) : (pbc.lattice.b,pbc.lattice.a)
        if !isinf(R)
            ϕ = pbc.shifts[ind]*K*R
        else
            ϕ = zeros(Float64,size(pbc.shifts[ind]))
        end
        return exp.(-1im*ϕ)
    end

    Base.show(io::IO,::FloquetBC{Union{},Union{},false}) = print("FloquetBC()")
    function Base.show(io::IO,PBC::FloquetBC{DIM,SIDE,FIN}) where {FIN,DIM,SIDE}
        print("FloquetBC{$DIM,$SIDE}")
        try
            print(",v1=",fmt("3.1f",PBC.lattice.a),"∠",fmt("3.1f",PBC.lattice.α*180/π))
            print(",v2=",fmt("3.1f",PBC.lattice.b),"∠",fmt("3.1f",PBC.lattice.β*180/π))
        catch
            print(" w/ unintialized lattice")
        end
    end
end


"""
    MatchedBC{Polar}(dim, side, N)
"""
struct MatchedBC{DIM,SIDE,CS,TBC1,TBC2,FIN} <: AbstractBC{DIM,SIDE}
    outgoing_qns::Array{Int,1}
    incoming_qns::Array{Int,1}
    N::NTuple{2,Int}
    dx::NTuple{2,Float64}
    xmin::NTuple{2,Float64}
    xmax::NTuple{2,Float64}
    I1::Array{Array{Int,1},1}
    J1::Array{Array{Int,1},1}
    I2::Array{Array{Int,1},1}
    J2::Array{Array{Int,1},1}
    weights::Array{Array{ComplexF64,1},1}
    QNs::Array{Int,1}
    direction::Array{Int,1}

    MatchedBC() = new{Union{},Union{},Union{},Union{},Union{},false}()
    MatchedBC{DIM,SIDE}(outgoing_qns, incoming_qns; kwargs...) where {DIM,SIDE} = new{DIM,SIDE,Union{},Union{},Union{},false}(outgoing_qns, incoming_qns)

    # MatchedBC for 2-dim free space polar
    function MatchedBC{1,2,Polar}(outgoing_qns, incoming_qns, N; kwargs...)
        common = intersect(outgoing_qns,incoming_qns)
        @assert isempty(common) "outgoing and incoming sets can share no common channel, but their intersection is $common"

        @assert haskey(kwargs,:dx) "MatchedBC requires kwarg :dx"
        dx = kwargs[:dx]
        xmin = haskey(kwargs,:xmin) ? kwargs[:xmin] : (0.0,0.0)
        xmax = (xmin[1]+(N[1]-1/2)*dx[1],xmin[2]+(N[2]-1/2)*dx[2])

        dr, dθ = dx[1], dx[2]

        q_out = length(outgoing_qns)
        q_in = length(incoming_qns)

        weights = Array{ComplexF64}(undef,(q_out+q_in)*N[2]^2)
        M  = Array{Int}(undef,q_out+q_in)
        direction = Array{Int}(undef,q_out+q_in)
        I1 = Array{Int}(undef,(q_out+q_in)*N[2]^2)
        I2 = Array{Int}(undef,(q_out+q_in)*N[2]^2)
        J1 = Array{Int}(undef,(q_out+q_in)*N[2]^2)
        J2 = Array{Int}(undef,(q_out+q_in)*N[2]^2)
        M_inds = Array{Int}(undef,(q_out+q_in)*N[2]^2)
        M_subinds = Array{Array{Int,1},1}(undef,(q_out+q_in))
        idx = 1
        midx = 1
        for m ∈ outgoing_qns
            for i ∈ 1:N[2], j ∈ 1:N[2]
                I1[idx] = N[1]
                J1[idx] = i
                I2[idx] = N[1]
                J2[idx] = j
                # weights[idx] = dθ*exp(complex(0,m*(i-j)*dθ))/2π
                weights[idx] = dθ*exp(complex(0,-m*(i-j)*dθ))/2π
                M_inds[idx] = midx
                idx += 1
            end
            M[midx] = m
            M_subinds[midx] = (midx-1)*N[2]^2 .+ (1:N[2]^2)
            direction[midx] = 1
            midx +=1
        end
        qidx = midx-1
        for m ∈ incoming_qns
        for i ∈ 1:N[2], j ∈ 1:N[2]
                I1[idx] = N[1]
                J1[idx] = i
                I2[idx] = N[1]
                J2[idx] = j
                # weights[idx] = dθ*exp(complex(0,m*(i-j)*dθ))/2π
                weights[idx] = dθ*exp(complex(0,-m*(i-j)*dθ))/2π
                M_inds[idx] = midx
                idx += 1
            end
            M[midx] = m
            M_subinds[midx] = qidx*N[2]^2 .+ (midx-qidx-1)*N[2]^2 .+ (1:N[2]^2)
            direction[midx] = 2
            midx += 1
        end
        i1 = Array{Array{Int,1},1}(undef,length(M))
        i2 = Array{Array{Int,1},1}(undef,length(M))
        j1 = Array{Array{Int,1},1}(undef,length(M))
        j2 = Array{Array{Int,1},1}(undef,length(M))
        v = Array{Array{ComplexF64,1},1}(undef,length(M))
        for i ∈ eachindex(M)
            i1[i] = I1[M_subinds[i]]
            i2[i] = I2[M_subinds[i]]
            j1[i] = J1[M_subinds[i]]
            j2[i] = J2[M_subinds[i]]
            v[i] = weights[M_subinds[i]]
        end
        return new{1,2,Polar,true,Union{},Union{}}(outgoing_qns, incoming_qns, Tuple(N), Tuple(dx), Tuple(xmin), xmax, i1, j1, i2, j2, v, M, direction)
    end
    MatchedBC{1,2,Polar,BC1,BC2}(args...; kwargs...) where {BC1,BC2} = MatchedBC{1,2,Polar}(args...; kwargs...)

    # MatchedBC for Cartesian waveguide system (not index guided)
    function MatchedBC{DIM,SIDE,Cartesian,BC1,BC2}(outgoing_qns, incoming_qns, N; kwargs...) where {DIM,SIDE,BC1,BC2}
        common = intersect(outgoing_qns,incoming_qns)
        @assert isempty(common) "outgoing and incoming sets can share no common channel, but their intersection is $common"
        @assert haskey(kwargs,:dx) "MatchedBC requires kwarg :dx"
        dim_perp = mod1(DIM+1,2)
        dx = kwargs[:dx]
        xmin = haskey(kwargs,:xmin) ? kwargs[:xmin] : nothing#(0.0,0.0)
        xmax = (xmin[1]+N[1]*dx[1],xmin[2]+N[2]*dx[2])
        DX = dx[dim_perp]
        L = (xmax.-xmin)[dim_perp]

        q_out = length(outgoing_qns)
        q_in = length(incoming_qns)

        M  = Array{Int}(undef,q_out+q_in)
        direction = Array{Int}(undef,q_out+q_in)
        weights = Array{ComplexF64}(undef,(q_out+q_in)*N[dim_perp]^2)
        I1 = Array{Int}(undef,(q_out+q_in)*N[dim_perp]^2)
        I2 = Array{Int}(undef,(q_out+q_in)*N[dim_perp]^2)
        J1 = Array{Int}(undef,(q_out+q_in)*N[dim_perp]^2)
        J2 = Array{Int}(undef,(q_out+q_in)*N[dim_perp]^2)
        M_inds = Array{Int}(undef,(q_out+q_in)*N[dim_perp]^2)
        M_subinds = Array{Array{Int,1},1}(undef,(q_out+q_in))
        idx = 1
        midx = 1
        for m ∈ outgoing_qns
            am, bm, km = abkm(m,L,BC1,BC2)
            for i ∈ 1:N[dim_perp], j ∈ 1:N[dim_perp]
                if SIDE==1
                    I1[idx] = 1
                    J1[idx] = i
                    I2[idx] = 1
                    J2[idx] = j
                else
                    I1[idx] = N[DIM]
                    J1[idx] = i
                    I2[idx] = N[DIM]
                    J2[idx] = j
                end
                ϕi = am*exp(1im*km*DX*(1/2+i-1))+bm*exp(-1im*km*DX*(1/2+i-1))
                ϕj = am*exp(1im*km*DX*(1/2+j-1))+bm*exp(-1im*km*DX*(1/2+j-1))
                weights[idx] = DX*ϕi*ϕj
                M_inds[idx] = midx
                idx += 1
            end
            M[midx] = m
            M_subinds[midx] = (midx-1)*N[dim_perp]^2 .+ (1:N[dim_perp]^2)
            direction[midx] = SIDE==1 ? 1 : 1
            midx +=1
        end
        qidx = midx-1
        for m ∈ incoming_qns
            am, bm, km = abkm(m,L,BC1,BC2)
            for i ∈ 1:N[dim_perp], j ∈ 1:N[dim_perp]
                if SIDE==1
                    I1[idx] = 1
                    J1[idx] = i
                    I2[idx] = 1
                    J2[idx] = j
                else
                    I1[idx] = N[DIM]
                    J1[idx] = i
                    I2[idx] = N[DIM]
                    J2[idx] = j
                end
                ϕi = am*exp(1im*km*DX*(1/2+i-1))+bm*exp(-1im*km*DX*(1/2+i-1))
                ϕj = am*exp(1im*km*DX*(1/2+j-1))+bm*exp(-1im*km*DX*(1/2+j-1))
                weights[idx] = DX*ϕi*ϕj
                M_inds[idx] = midx
                idx += 1
            end
            M[midx] = m
            M_subinds[midx] = qidx*N[dim_perp]^2 .+ (midx-qidx-1)*N[dim_perp]^2 .+ (1:N[dim_perp]^2)
            direction[midx] = SIDE==1 ? -1 : -1
            midx += 1
        end
        i1 = Array{Array{Int,1},1}(undef,length(M))
        i2 = Array{Array{Int,1},1}(undef,length(M))
        j1 = Array{Array{Int,1},1}(undef,length(M))
        j2 = Array{Array{Int,1},1}(undef,length(M))
        v = Array{Array{ComplexF64,1},1}(undef,length(M))
        for i ∈ eachindex(M)
            i1[i] = I1[M_subinds[i]]
            i2[i] = I2[M_subinds[i]]
            j1[i] = J1[M_subinds[i]]
            j2[i] = J2[M_subinds[i]]
            v[i] = weights[M_subinds[i]]
        end
        return new{DIM,SIDE,Cartesian,BC1,BC2,true}(outgoing_qns, incoming_qns, Tuple(N), Tuple(dx), Tuple(xmin), xmax, i1, j1, i2, j2, v, M, direction)
    end
    # MatchedBC{DIM,SIDE,Cartesian,BC1,BC2}(args...; kwargs...) where {DIM,SIDE,BC1,BC2} = MatchedBC{DIM,SIDE,Cartesian}(args...; kwargs...)

    # convenience constructors for MatchedBC, starting from partially initialized
    # (mbc::MatchedBC{DIM,SIDE,_1,_2,_3,false})(N, coordinate_system::Type{CS}, bc1::Type{BC1}, bc2::Type{BC2}, args...; kwargs...) where {DIM,SIDE,_1,_2,_3,CS,BC1,BC2} = MatchedBC{DIM,SIDE,CS,BC1,BC2}(mbc.outgoing_qns,mbc.incoming_qns,N; kwargs...)
    (mbc::MatchedBC{DIM,SIDE,CS,BC1,BC2,false})(N) where {DIM,SIDE,CS,BC1,BC2}= MatchedBC{DIM,SIDE,CS,BC1,BC2}(mbc.outgoing_qns,mbbc.incoming_qns,N)

    # longitudinal propagator for free polar
    function (MBC::MatchedBC{1,2,Polar})(k,args...)
        R = MBC.xmax[1]+MBC.dx[1]/2

        Hm = besselhx.(MBC.QNs,MBC.direction,k*R)

        Hm₋ = besselhx.(MBC.QNs.-1,MBC.direction,k*R)
        Hm₊ = besselhx.(MBC.QNs.+1,MBC.direction,k*R)
        Hm′ = (Hm₋-Hm₊)/2

        αm = Hm′./Hm
        kdR = k*MBC.dx[1]
        num = 1 .+ αm*kdR/2
        den = 1 .- αm*kdR/2
        return num./den
    end
    # longitudinal propagator for waveguide cartesian
    function (MBC::MatchedBC{DIM,SIDE,Cartesian,BC1,BC2})(k,args...) where {DIM,SIDE,BC1,BC2}
        km = Array{Float64}(undef,length(MBC.QNs))
        for i ∈ eachindex(MBC.QNs)
            _, _, km[i] = abkm.(MBC.QNs[i],(MBC.xmax.-MBC.xmin)[mod1(DIM+1,2)],BC1,BC2)
        end
        β = @. MBC.direction*sqrt(k^2-km^2 + 0.0im)
        βdx = β*MBC.dx[DIM]#[mod1(DIM+1,2)]
        num = 1 .+ 1im*βdx/2
        den = 1 .- 1im*βdx/2
        return num./den
    end

    # pretty printing for MatchedBC
    Base.show(io::IO,::MatchedBC{Union{},Union{}}) = print("MatchedBC()")
    Base.show(io::IO,::MatchedBC{DIM,SIDE}) where {DIM,SIDE} = print("MatchedBC{$DIM,$SIDE}")
    Base.show(io::IO,::MatchedBC{DIM,SIDE,CS}) where {DIM,SIDE,CS<:CoordinateSystem} = print("MatchedBC{$DIM,$SIDE,$CS}")
end
# coefficients and transverse momentum for Cartesian guided mode system
function abkm(m::Int,L::Real,bc1::Type{BC1},bc2::Type{BC2}) where {BC1,BC2}
    if BC1<:DirichletBC && BC2<:DirichletBC
        am =  sqrt(2/L)/2im
        bm = -sqrt(2/L)/2im
        km = m*π/L
    elseif BC1<:DirichletBC && BC2<:NeumannBC
        throw(ErrorException("haven't done Dirichlet, Neumann pairing yet"))
        am =  sqrt(2/L)/2im
        bm = -sqrt(2/L)/2im
        km = m*π/L
    elseif BC1<:NeumannBC   && BC2<:DirichletBC
        throw(ErrorException("haven't done Dirichlet, Neumann pairing yet"))
        am =  sqrt(2/L)/2im
        bm = -sqrt(2/L)/2im
        km = m*π/L
    elseif BC1<:DirichletBC && BC2<:NeumannBC
        throw(ErrorException("haven't done Dirichlet, Neumann pairing yet"))
        am =  sqrt(2/L)/2im
        bm = -sqrt(2/L)/2im
        km = m*π/L
    elseif BC1<:NeumannBC   && BC2<:NeumannBC
        am = sqrt(2/L)/2
        bm = sqrt(2/L)/2
        km = m*π/L
    elseif BC1<:FloquetBC  && BC2<:FloquetBC
        am = m<0 ?  sqrt(2/L)/2im : sqrt(2/L)/2
        bm = m<0 ? -sqrt(2/L)/2im : sqrt(2/L)/2
        km = abs(m)*2π/L
    else
        throw(Exception("unrecognized transverse boundary condition combination $BC1 & $BC2"))
    end
    return am, bm, km
end

#################### DIM/SIDE UTILITIES #####################
get_bc_type(::DirichletBC) = Float64
get_bc_type(::NeumannBC) = Float64
get_bc_type(::RobinBC{DIM,SIDE,FIN,TW}) where {FIN,DIM,SIDE,TW<:Array} = eltype(TW)
get_bc_type(::RobinBC{DIM,SIDE,FIN,TW}) where {FIN,DIM,SIDE,TW} = TW
get_bc_type(::FloquetBC) = ComplexF64
get_bc_type(::MatchedBC) = ComplexF64

get_dim_side(x::Union{Tuple,Array}) = map(get_dim_side,x)
get_dim_side(x) = nothing
get_dim_side(::BC) where BC<:AbstractBC{DIM,SIDE} where {DIM,SIDE} = Val{DIM}(),Val{SIDE}()
get_dim_side(::BL) where BL<:AbstractBL{DIM,SIDE} where {DIM,SIDE} = Val{DIM}(),Val{SIDE}()

get_dim(x::Union{Tuple,Array}) = map(get_dim,x)
get_dim(x) = nothing
get_dim(BCL::Union{AbstractBC,AbstractBL}) = get_dim_side(BCL)[1]

get_side(x::Union{Tuple,Array}) = map(get_side,x)
get_side(x) = nothing
get_side(BCL::Union{AbstractBC,AbstractBL}) = get_dim_side(BCL)[2]

apply_args(x::Tuple, args...;kwargs...) = map(x->apply_args(x, args...;kwargs...),x)
# this next one is to help with type stability in the recursion
apply_args(x::Tuple{Tuple,Tuple}, args...;kwargs...) = map(x->apply_args(x, args...;kwargs...),x)
apply_args(x, args...;kwargs...) = nothing
apply_args(x::AbstractLocalBC, args...; N, kwargs...) =
    x(N; kwargs...)
apply_args(x::FloquetBC{DIM,SIDE,false}, args...; lattice, N, kwargs...) where {DIM,SIDE} =
    x(lattice, N; kwargs...)
apply_args(x::FloquetBC{DIM,SIDE,true}, args...; lattice, N, kwargs...) where {DIM,SIDE}= x
apply_args(x::MatchedBC{DIM,SIDE,CS1,BC1,BC2,false},
    coordinate_system::CS2, bc1::BC3, bc2::BC4, args... ;
    N, kwargs...) where {DIM,SIDE,CS1,CS2<:CoordinateSystem,BC1,BC2,BC3,BC4} =
        MatchedBC{DIM,SIDE,CS2,BC3,BC4}(x.outgoing_qns, x.incoming_qns, N; kwargs...)
apply_args(x::MatchedBC{DIM,SIDE,CS1,BC1,BC2,true}, args...; N, kwargs...) where {DIM,SIDE,CS1,BC1,BC2} = x
# apply_args(x::AbstractBL, args...; Δ, kwargs...) = x(Δ)
apply_args(x::noBL, args...; kwargs...) = x
apply_args(x::AbstractBL, args...; depth::Number=x.depth, Δ, bounds=x.bounds, kwargs...) = x(depth,Δ,bounds)

reorder_side(bcs::Tuple) = reorder_side(bcs...)
reorder_side(bc,bcs...) = typeof(get_side(bc))<:Val{1} ? (bc,reorder_side(bcs...)...) : (reorder_side(bcs...)...,bc)
reorder_side(bc) = (bc,)

reorder_dim(bcs::Tuple) = reorder_dim(bcs...)
reorder_dim(bc,bcs...) = typeof(get_dim(bc))<:Val{1} ? (bc,reorder_dim(bcs...)...) : (reorder_dim(bcs...)...,bc)
reorder_dim(bc) = (bc,)

reorder(bcs...) = reorder_dim(reorder_dim(reorder_side(bcs...)))



"""
    isDirichlet(bc)
"""
isDirichlet(bc::DirichletBC) = true
isDirichlet(bc) = false


"""
    isNeumann(bc)
"""
isNeumann(bc::NeumannBC) = true
isNeumann(bc) = false


"""
    isMatched(bc)
"""
isMatched(bc::MatchedBC) = true
isMatched(bc) = false


"""
    isFloquet(bc)
"""
isFloquet(bc::FloquetBC) = true
isFloquet(bc) = false


"""
    isPML(bl)
"""
isPML(bl::PML) = true
isPML(bl) = false


"""
    iscPML(bl)
"""
iscPML(bl::cPML) = true
iscPML(bl) = false


"""
    isnoBL(bl)
"""

isnoBL(bl::noBL) = true
isnoBL(bl) = false




# #################### EZBC/EZBL UTILITY ###########################
# """
#     bcs = ezbc(bc1, bc2, N=[0,0])
#     bcs = ezbc(bcs, N=[0,0])
#
# wrapper for easy construction of boundary conditions, `bc1` for dimension 1, `bc2` for dimension 2.
#
# `bc1`, `bc2` can by specified as symbols or arrays of symbols. If symbol, it applies to both boundaries in that dimension.
#
# See also: [`RobinBoundaryCondition`](@ref), [`PeriodicBoundaryCondition`](@ref)
# """
# ezbc(bc1::Symbol,bc2; kwargs...) = ezbc([bc1, bc1],bc2; kwargs...)
# ezbc(bc1,bc2::Symbol; kwargs...) =  ezbc(bc1, [bc2,bc2]; kwargs...)
# ezbc(bc1::Symbol, bc2::Symbol; kwargs...) = ezbc([bc1, bc1], [bc2, bc2]; kwargs...)
# ezbc(bc::Symbol; kwargs...) = ezbc(bc,bc; kwargs...)
# function ezbc(bc1,bc2; kwargs...)
#
#     dirichlet = [:d,:D,:dirichlet,:Dirichlet]
#     neumann = [:n,:N,:neumann,:Neumann]
#     robin = [:r,:R,:robin,:Robin]
#     periodic = [:p,:P,:periodic,:Periodic,:floquet,:Floquet]
#     matched = [:m,:M,:matched,:Matched]
#     valid_names = vcat(dirichlet,neumann,robin,periodic,matched)
#
#     @assert bc1[1] ∈ valid_names "unrecognized boundary condition for dim=1, side=1"
#     @assert bc1[2] ∈ valid_names "unrecognized boundary condition for dim=1, side=2"
#     @assert bc2[1] ∈ valid_names "unrecognized boundary condition for dim=2, side=1"
#     @assert bc2[2] ∈ valid_names "unrecognized boundary condition for dim=2, side=2"
#
#     bc1[1] ∈ dirichlet ? BC11 = DirichletBC{1,1}() : nothing
#     bc1[1] ∈ neumann ? BC11 = NeumannBC{1,1}() : nothing
#     bc1[1] ∈ matched ? BC11 = MatchedBC{1,1}(kwargs[:outgoing_qns],kwargs[:incoming_qns]) : nothing
#     bc1[1] ∈ periodic ? BC11 = FloquetBC{1,1}(kwargs[:lattice]) : nothing
#     bc1[1] ∈ robin ? BC11 = RobinBC{1,1}(kwargs[:α],kwargs[:β]) : nothing
#
#     bc1[2] ∈ dirichlet ? BC12 = DirichletBC{1,2}() : nothing
#     bc1[2] ∈ neumann ? BC12 = NeumannBC{1,2}() : nothing
#     bc1[2] ∈ matched ? BC12 = MatchedBC{1,2}(kwargs[:outgoing_qns],kwargs[:incoming_qns]) : nothing
#     bc1[2] ∈ periodic ? BC12 = FloquetBC{1,2}(kwargs[:lattice]) : nothing
#     bc1[2] ∈ robin ? BC12 = RobinBC{1,2}(kwargs[:α],kwargs[:β]) : nothing
#
#     bc2[1] ∈ dirichlet ? BC21 = DirichletBC{2,1}() : nothing
#     bc2[1] ∈ neumann ? BC21 = NeumannBC{2,1}() : nothing
#     bc2[1] ∈ matched ? BC21 = MatchedBC{2,1}(kwargs[:outgoing_qns],kwargs[:incoming_qns]) : nothing
#     bc2[1] ∈ periodic ? BC21 = FloquetBC{2,1}(kwargs[:lattice]) : nothing
#     bc2[1] ∈ robin ? BC21 = RobinBC{2,1}(kwargs[:α],kwargs[:β]) : nothing
#
#     bc2[2] ∈ dirichlet ? BC22 = DirichletBC{2,2}() : nothing
#     bc2[2] ∈ neumann ? BC22 = NeumannBC{2,2}() : nothing
#     bc2[2] ∈ matched ? BC22 = MatchedBC{2,2}(kwargs[:outgoing_qns],kwargs[:incoming_qns]) : nothing
#     bc2[2] ∈ periodic ? BC22 = FloquetBC{2,2}(kwargs[:lattice]) : nothing
#     bc2[2] ∈ robin ? BC22 = RobinBC{2,2}(kwargs[:α],kwargs[:β]) : nothing
#
#     return ((BC11,BC12),(BC21,BC22))
# end
#
#
# ezbl(bl1::Symbol,bl2; kwargs...) = ezbl([bl1, bl1],bl2; kwargs...)
# ezbl(bl1,bl2::Symbol; kwargs...) = ezbl(bl1, [bl2,bl2]; kwargs...)
# ezbl(bl1::Symbol, bl2::Symbol; kwargs...) = ezbl([bl1, bl1], [bl2, bl2]; kwargs...)
# ezbl(bl::Symbol; kwargs...) = ezbl(bl,bl; kwargs...)
# function ezbl(bl1,bl2; kwargs...)
#
#     pml = [:p,:P,:pml,:PML,:pml_out,:PML_OUT,:PML_out]
#     cpml = [:c,:C,:cp,:CP,:cP,:cpml,:cPML,:conj_pml,:conj_PML,:pml_in,:PML_IN,:PML_in]
#     nobl = [:n,:none,:N,:None,:NONE,:no_bl]
#     valid_names = vcat(pml,cpml,nobl)
#     @assert bl1[1] ∈ valid_names "unrecognized boundary layer for dim=1, side=1"
#     @assert bl1[2] ∈ valid_names "unrecognized boundary layer for dim=1, side=2"
#     @assert bl2[1] ∈ valid_names "unrecognized boundary layer for dim=2, side=1"
#     @assert bl2[2] ∈ valid_names "unrecognized boundary layer for dim=2, side=2"
#
#     bl1[1] ∈ pml ? BL11 = PML{1,1}() : nothing
#     bl1[1] ∈ cpml ? BL11 = cPML{1,1}() : nothing
#     bl1[1] ∈ nobl ? BL11 = noBL{1,1}() : nothing
#
#     bl1[2] ∈ pml ? BL12 = PML{1,2}() : nothing
#     bl1[2] ∈ cpml ? BL12 = cPML{1,2}() : nothing
#     bl1[2] ∈ nobl ? BL12 = noBL{1,2}() : nothing
#
#     bl2[1] ∈ pml ? BL21 = PML{2,1}() : nothing
#     bl2[1] ∈ cpml ? BL21 = cPML{2,1}() : nothing
#     bl2[1] ∈ nobl ? BL21 = noBL{2,1}() : nothing
#
#     bl2[2] ∈ pml ? BL22 = PML{2,2}() : nothing
#     bl2[2] ∈ cpml ? BL22 = cPML{2,2}() : nothing
#     bl2[2] ∈ nobl ? BL22 = noBL{2,2}() : nothing
#
#     return ((BL11,BL12),(BL21,BL22))
# end

end # module
