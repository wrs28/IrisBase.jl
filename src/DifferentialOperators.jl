# TODO: check effect of polar coordinates on ihomogeneous systems.
# TODO: loop over size of M to apply BC's deeper in the interior in order to do photonic crystal waveguides as open boundaries
# TODO: fix gradient

module DifferentialOperators

using ..Bravais,
..BoundaryConditions,
..CoordinateSystems,
LinearAlgebra,
SparseArrays

import ..BoundaryConditions: get_dim_side, get_dim, get_side, get_bc_type, reorder, apply_args

export laplacian


struct laplacian{DIM,POS} end
"""
    laplacian(N, coordinate_system, Δ, bcs[, bls, k, ka, kb, ind]; kwargs...)

    laplacian{DIM}(N, coordinate_system, Δ, bcs[, bls, k, ka, kb, ind]; kwargs...)
"""
function laplacian(N, coordinate_system::Type{TCS}, Δ, bcs, bls=ezbl(:n), args...; kwargs...) where TCS<:CoordinateSystem
    OD1 = define_op{2,1,TCS}(N,Δ,bcs,bls; kwargs...)
    OD2 = define_op{2,2,TCS}(N,Δ,bcs,bls; kwargs...)
    ∇₁² = diff_op(OD1,args...)
    ∇₂² = diff_op(OD2,args...)
    return sparse(∇₁², ∇₂²), (OD1,OD2)
end
function laplacian{DIM}(N, coordinate_system::Type{TCS}, Δ, bcs, bls=ezbl(:n), args...; kwargs...) where {DIM,TCS<:CoordinateSystem}
    N = DIM==1 ? (N[1],1) : (1,N[2])
    OD = define_op{2,DIM,TCS}(N,Δ,bcs,bls; kwargs...)
    ∇² = diff_op(OD,args...)
    return sparse(∇²), (OD,)
end

##############################################################################################
abstract type AbstractPolarity end
struct Backward<:AbstractPolarity end
struct Forward<:AbstractPolarity end
struct Central<:AbstractPolarity end


struct OperatorDefinition{ORD,DIM,TCS,TBC11,TBC12,TBC21,TBC22,TBL11,TBL12,TBL21,TBL22,TH1,TH2,TF,POL}
    N::NTuple{2,Int}
    Δ::NTuple{2,NTuple{2,Float64}}
    dx::NTuple{2,Float64}
    x::NTuple{2,LinRange{Float64}}
    bcs::Tuple{Tuple{TBC11,TBC12},Tuple{TBC21,TBC22}}
    bls::Tuple{Tuple{TBL11,TBL12},Tuple{TBL21,TBL22}}
    h::Tuple{Array{TH1,1},Array{TH2,1}}
    f::NTuple{2,Array{TF,1}}
    xd::NTuple{2,LinRange{Float64}}
    hd::Tuple{Array{TH1,1},Array{TH2,1}}
    fd::NTuple{2,Array{TF,1}}

    OperatorDefinition{ORD,DIM,CS}(N,Δ,bcs; kwargs...) where {ORD,DIM,CS} = OperatorDefinition{ORD,DIM,CS}(N,_oc_Δ(Δ),_oc_bcs(bcs),((noBL{1,1}(),noBL{1,2}()),(noBL{2,1}(),noBL{2,2}())); kwargs...)
    OperatorDefinition{ORD,DIM,CS}(N,Δ,bcs,bls; kwargs...) where {ORD,DIM,CS} = OperatorDefinition{ORD,DIM,CS}(N,_oc_Δ(Δ),_oc_bcs(bcs),_oc_bls(bls); kwargs...)
    function OperatorDefinition{ORD,DIM,CS}(N,Δ::Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}},bcs::Tuple{Tuple,Tuple},bls::Tuple{Tuple,Tuple}; kwargs...) where {ORD,DIM,CS,TD}

        bc1_bool = typeof(bcs[1][1])<:PeriodicBC && typeof(bcs[1][2])<:PeriodicBC
        bc2_bool = typeof(bcs[2][1])<:PeriodicBC && typeof(bcs[2][2])<:PeriodicBC
        if (bc1_bool || bc2_bool)
            lattice = haskey(kwargs,:lattice) ? kwargs[:lattice] : BravaisLattice()
        else
            lattice = BravaisLattice()
        end
        a, b, α, β = lattice.a, lattice.b, lattice.α, lattice.β

        dx1 = !isinf(a) ? a*sin(β-α)/sin(β)/N[1] : (Δ[1][2]-Δ[1][1])/N[1]
        dx2 = !isinf(b) ?  b*sin(β)/N[2] : (Δ[2][2]-Δ[2][1])/N[2]

        Δ = !isinf(a) ? ((Δ[1][1],Δ[1][1]+a*sin(β-α)/sin(β)),Δ[2]) : Δ
        Δ = !isinf(b) ? (Δ[1],(Δ[2][1],Δ[2][1]+b*sin(β))) : Δ

        xmin = float.((Δ[1][1],Δ[2][1]))
        dx = float.((dx1,dx2))

        x1 = LinRange(Δ[1][1]+dx[1]/2,Δ[1][2]-dx[1]/2,N[1])
        x2 = LinRange(Δ[2][1]+dx[2]/2,Δ[2][2]-dx[2]/2,N[2])
        x = (x1,x2)
        x1 = LinRange(Δ[1][1],Δ[1][2],N[1]+1)
        x2 = LinRange(Δ[2][1],Δ[2][2],N[2]+1)
        xd = (x1,x2)

        bcs = apply_args(bcs, CS ; N=N, dx=dx, xmin=xmin, lattice=lattice, kwargs...)
        depth = haskey(kwargs,:depth) ? kwargs[:depth] : 0
        bls = apply_args(bls; Δ=Δ, depth=depth)
        h1 = typeof(bls[1])<:Tuple{noBL,noBL} ? zeros(Int,size(x[1])) : 1im*(bls[1][1].(x[1]) + bls[1][2].(x[1]))
        h2 = typeof(bls[2])<:Tuple{noBL,noBL} ? zeros(Int,size(x[2])) : 1im*(bls[2][1].(x[2]) + bls[2][2].(x[2]))
        h = (h1,h2)
        h1 = typeof(bls[1])<:Tuple{noBL,noBL} ? zeros(Int,size(xd[1])) : 1im*(bls[1][1].(xd[1]) + bls[1][2].(xd[1]))
        h2 = typeof(bls[2])<:Tuple{noBL,noBL} ? zeros(Int,size(xd[2])) : 1im*(bls[2][1].(xd[2]) + bls[2][2].(xd[2]))
        hd = (h1,h2)
        f1 = typeof(bls[1])<:Tuple{noBL,noBL} ? zeros(Int,size(x[1])) : 1im*(-reverse(cumsum(bls[1][1].(x[1][end:-1:1]))) + cumsum(bls[1][2].(x[1])))*dx[1]
        f2 = zeros(Int,size(x[2]))
        f = (f1,f2)
        fd1 = typeof(bls[1])<:Tuple{noBL,noBL} ? zeros(Int,size(xd[1])) : 1im*(-reverse(cumsum(bls[1][1].(xd[1][end:-1:1]))) + cumsum(bls[1][2].(xd[1])))*dx[1]
        fd2 = zeros(Int,size(xd[2]))
        fd = (fd1,fd2)

        pol = haskey(kwargs,:polarity) ? typeof(kwargs[:polarity]) : (ORD==1 ? Central : Union{})

        tbc11 = typeof(bcs[1][1])
        tbc12 = typeof(bcs[1][2])
        tbc21 = typeof(bcs[2][1])
        tbc22 = typeof(bcs[2][2])

        tbl11 = typeof(bls[1][1])
        tbl12 = typeof(bls[1][2])
        tbl21 = typeof(bls[2][1])
        tbl22 = typeof(bls[2][2])

        return new{ORD,DIM,CS,tbc11,tbc12,tbc21,tbc22,tbl11,tbl12,tbl21,tbl22,eltype(h1),eltype(h2),eltype(f1),pol}(Tuple(N),Δ,dx,x,bcs,bls,h,f,xd,hd,fd)
    end

    Base.show(io::IO,OD::OperatorDefinition{ORD,DIM,TCS}) where {ORD,DIM,TCS} = print("OperatorDefinition{$ORD,$DIM,$TCS}")
end
_oc_Δ(Δ::Array{T,2}) where T = _oc_Δ(((Δ[1,1],Δ[2,1]),(Δ[1,2],Δ[2,2])))
_oc_Δ(Δ::Array{T,1}) where T = _oc_Δ(((Δ[1],Δ[2]),(Δ[3],Δ[4])))
_oc_Δ(Δ::Tuple) = _oc_Δ(((Δ[1],Δ[2]),(Δ[3],Δ[4])))
_oc_Δ(Δ::Tuple{Tuple,Tuple}) = ((float(Δ[1][1]),float(Δ[1][2])),(float(Δ[2][1]),float(Δ[2][2])))
_oc_Δ(Δ) = Δ

function _oc_bcs(bcs)
    bcs = reorder(bcs)
    return ((bcs[1],bcs[2]),(bcs[3],bcs[4]))
end
_oc_bcs(bcs::Tuple{Tuple,Tuple}) = bcs

_oc_bls(bls::Tuple) = _oc_bls(bls...)
function _oc_bls(bl...)
    bl = length(bl)==0 ? reorder(bl...,noBL{1,1}(),) : reorder(bl)
    bl = typeof(bl[1])<:AbstractBL{1,1} ? (bl...,) : (noBL{1,1}(),bl...)

    bl = length(bl)==1 ? reorder(bl...,noBL{1,2}(),) : reorder(bl)
    bl = typeof(bl[2])<:AbstractBL{1,2} ? (bl...,) : (noBL{1,2}(),bl...)

    bl = length(bl)==2 ? reorder(bl...,noBL{2,1}(),) : reorder(bl)
    bl = typeof(bl[3])<:AbstractBL{2,1} ? (bl...,) : (noBL{2,1}(),bl...)

    bl = length(bl)==3 ? reorder(bl...,noBL{2,2}(),) : reorder(bl)
    bl = typeof(bl[4])<:AbstractBL{2,2} ? (bl...,) : (noBL{2,2}(),bl...)
    return ((bl[1],bl[2]),(bl[3],bl[4]))
end
_oc_bls(bls::Tuple{Tuple,Tuple}) = bls
struct define_op{ORD,DIM,CS} end

"""
    define_op{ORD,DIM,CS}(args...; kwargs...)
"""
define_op{ORD,DIM,CS}(args...; kwargs...) where {ORD,DIM,CS} = OperatorDefinition{ORD,DIM,CS}(args...;kwargs...)

##################################################################

struct RowColVal{T}
    N::NTuple{2,Int}
    str::NTuple{2,Int}
    rows::Array{Int,1}
    cols::Array{Int,1}
    vals::Array{T,1}

    RowColVal(N,rows,cols,vals) = new{eltype(vals)}(N,Base.size_to_strides(1,N...),rows,cols,vals)
    RowColVal(N) = new{Union{}}(Tuple(N),Base.size_to_strides(1,N...))
    RowColVal(N,vec) = new{eltype(vec)}(N,Base.size_to_strides(1,N...),1:prod(N),1:prod(N),vec)
    function RowColVal(N,vec,dim)
        vec = dim==1 ? repeat(vec;outer=N[2]) : repeat(vec;inner=N[1])
        new{eltype(vec)}(N,Base.size_to_strides(1,N...),1:prod(N),1:prod(N),vec)
    end
    function (rcv::RowColVal{DIM})(i1, j1, i2, j2, vals) where DIM
        rows = 1 .+ rcv.str[1]*(i1.-1) + rcv.str[2]*(j1.-1)
        cols = 1 .+ rcv.str[1]*(i2.-1) + rcv.str[2]*(j2.-1)
        return new{eltype(vals)}(rcv.N,rcv.str,rows,cols,vals)
    end
    RowColVal(N,S::SparseMatrixCSC) = RowColVal(N,findnz(S)...)

    SparseArrays.sparse(rcv::RowColVal) = dropzeros!(sparse(rcv.rows,rcv.cols,rcv.vals,prod(rcv.N),prod(rcv.N),+))

    Base.promote_rule(::Type{RowColVal{T}}, ::Type{RowColVal{U}}) where {T,U} = RowColVal{promote_rule(T,U)}
    Base.convert(::Type{RowColVal{T}},rcv::RowColVal{U}) where {T,U} = RowColVal(rcv.N,rcv.rows,rcv.cols,convert(Array{T,1},rcv.vals))
    Base.convert(::Type{RowColVal{T}},rcv::RowColVal{T}) where T = rcv
end

##################################################################################

struct OperatorConstructor{ORD,DIM,CS,TBL1,TBL2,TBK,TB1,TB2,TOD}
    N::NTuple{2,Int}
    Δ::NTuple{2,NTuple{2,Float64}}
    dx::NTuple{2,Float64}
    bl1::TBL1
    bl2::TBL2
    bulk::RowColVal{TBK}
    bnd1::Array{RowColVal{TB1},1}
    bnd2::Array{RowColVal{TB2},1}
    definition::TOD

    function OperatorConstructor(OD::OperatorDefinition{ORD,DIM,CS,_1,_2,_3,_4,TBL1,TBL2,TBL3,TBL4}) where {ORD,DIM,CS,_1,_2,_3,_4,TBL1,TBL2,TBL3,TBL4}
        bulk, bnd1, bnd2 = RowColVal(OD.N), [RowColVal(OD.N)], [RowColVal(OD.N)]
        DIM==1 ? (tbl1,tbl2)=(TBL1,TBL2) : (tbl1,tbl2)=(TBL3,TBL4)
        return new{ORD,DIM,CS,tbl1,tbl2,Union{},Union{},Union{},typeof(OD)}(OD.N,OD.Δ,OD.dx,OD.bls[DIM][1],OD.bls[DIM][2],bulk,bnd1,bnd2,OD)
    end

    function (OC::OperatorConstructor{2,DIM,CS,_1,_2,Q,_3,_4,TZ})(args::Vararg{Number,M}) where {M,DIM,CS,_1,_2,_3,_4,Q,TZ}
        N,Δ,dx,def = OC.N,OC.Δ,OC.dx,OC.definition

        bulk_rows_on_site = 2:N[DIM]-1
        bulk_rows_nn = 1:N[DIM]-1
        bnd1_rows = 1
        bnd2_rows = N[DIM]
        rows = vcat(bulk_rows_on_site,bulk_rows_nn,bulk_rows_nn.+1,bnd1_rows,bnd2_rows)
        cols = vcat(bulk_rows_on_site,bulk_rows_nn.+1,bulk_rows_nn,bnd1_rows,bnd2_rows)
        vals = get_vals(OC,args...)

        if CS<:Polar && DIM==2
            rf = def.x[1] .+ def.f[1]/args[1]
            rf⁻² = 1 ./rf.^2
            S = spdiagm(0=>rf⁻²)
        else
            S = DIM==1 ? sparse(1.0*I,N[2],N[2]) : sparse(1.0*I,N[1],N[1])
        end
        rows,cols,vals = DIM==1 ? findnz(kron(S,sparse(rows,cols,vals,N[1],N[1],+))) : findnz(kron(sparse(rows,cols,vals,N[2],N[2],+),S))
        rcv = RowColVal(N,rows,cols,vals)
        return new{2,DIM,CS,_1,_2,eltype(vals),_3,_4,TZ}(OC.N,OC.Δ,OC.dx,OC.bl1,OC.bl2,rcv,OC.bnd1,OC.bnd2,OC.definition)
    end

    # different dims
    (OC::OperatorConstructor{ORD,DIM1})(BC::AbstractBC{DIM2},args...) where {ORD,DIM1,DIM2} = OC

    # same dims, side=1
    function (OC::OperatorConstructor{_1,DIM,CS,_3,_4,_5,_6,_7,_8})(BC::AbstractBC{DIM,1},args...) where {_1,DIM,CS,_3,_4,_5,_6,_7,_8}
        if typeof(BC)<:PeriodicBC
            @assert length(args)==4 "only $(length(args)) arguments passed. PeriodicBC necessitates 4 args: k, ka, kb, ind"
            # ind_flag = args[4]==0 ? false : true
            # ind_flag ? nothing : args[4]=1
            # ind = args[4]
            args[4]==0 ? args=(args[1:3]...,1) : nothing
            ind = args[4]
            i1,j1,i2,j2 = BC.I1[ind],BC.J1[ind],BC.I2[ind],BC.J2[ind]
        else
            i1,j1,i2,j2 = BC.I1,BC.J1,BC.I2,BC.J2
        end

        vals = get_vals(OC,BC,args...)
        rcv = RowColVal(OC.N)
        rcv = rcv.(i1, j1, i2, j2, vals)

        N = OC.N
        if CS<:Polar && DIM==2
            rf = OC.definition.x[1] .+ OC.definition.f[1]/args[1]
            rf⁻² = 1 ./rf.^2
            S = kron(sparse(I,N[2],N[2]),spdiagm(0=>rf⁻²))
            RCV = Array{RowColVal{eltype(S)},1}(undef,length(rcv))
        else
            S = sparse(I,prod(N),prod(N))
            RCV = rcv
        end
        for i ∈ eachindex(rcv)
            RCV[i] = RowColVal(N,S*sparse(rcv[i]))
        end
        return new{_1,DIM,CS,_3,_4,_5,eltype(vals[1]),_7,_8}(OC.N,OC.Δ,OC.dx,OC.bl1,OC.bl2,OC.bulk,RCV,OC.bnd2,OC.definition)
    end
    # same dims, side=2
    function (OC::OperatorConstructor{_1,DIM,CS,_3,_4,_5,_6,_7,_8})(BC::AbstractBC{DIM,2},args...) where {_1,DIM,CS,_3,_4,_5,_6,_7,_8}
        if typeof(BC)<:PeriodicBC
            @assert length(args)==4 "only $(length(args)) arguments passed. PeriodicBC necessitates 4 args: k, ka, kb, ind"
            args[4]==0 ? args=(args[1:3]...,1) : nothing
            ind = args[4]
            i1,j1,i2,j2 = BC.I1[ind],BC.J1[ind],BC.I2[ind],BC.J2[ind]
        else
            i1,j1,i2,j2 = BC.I1,BC.J1,BC.I2,BC.J2
        end

        vals = get_vals(OC,BC,args...)
        rcv = RowColVal(OC.N)
        rcv = rcv.(i1, j1, i2, j2, vals)

        N = OC.N
        if CS<:Polar && DIM==2
            rf = OC.definition.x[1] .+ OC.definition.f[1]/args[1]
            rf⁻² = 1 ./rf.^2
            S = kron(sparse(I,N[2],N[2]),spdiagm(0=>rf⁻²))
            RCV = Array{RowColVal{eltype(S)},1}(undef,length(rcv))
        else
            S = sparse(I,prod(N),prod(N))
            RCV = rcv
        end
        for i ∈ eachindex(rcv)
            RCV[i] = RowColVal(N,S*sparse(rcv[i]))
        end
        return new{_1,DIM,CS,_3,_4,_5,_6,eltype(vals[1]),_8}(OC.N,OC.Δ,OC.dx,OC.bl1,OC.bl2,OC.bulk,OC.bnd1,RCV,OC.definition)
    end

    Base.show(io::IO,::OperatorConstructor{ORD,DIM,_1,_2,_3,_4,Union{},Union{},Union{}}) where {ORD,DIM,_1,_2,_3,_4} = print("OC{$ORD,$DIM}: 0/3")
    Base.show(io::IO,::OperatorConstructor{ORD,DIM,_1,_2,_3,_4,TBK,Union{},Union{}}) where {ORD,DIM,_1,_2,_3,_4,TBK} = print("OC{$ORD,$DIM}: 1/3")
    Base.show(io::IO,::OperatorConstructor{ORD,DIM,_1,_2,_3,_4,Union{},TB1,Union{}}) where {ORD,DIM,_1,_2,_3,_4,TB1} = print("OC{$ORD,$DIM}: 1/3")
    Base.show(io::IO,::OperatorConstructor{ORD,DIM,_1,_2,_3,_4,Union{},Union{},TB2}) where {ORD,DIM,_1,_2,_3,_4,TB2} = print("OC{$ORD,$DIM}: 1/3")
    Base.show(io::IO,::OperatorConstructor{ORD,DIM,_1,_2,_3,_4,TBK,TB1,Union{}}) where {ORD,DIM,_1,_2,_3,_4,TBK,TB1} = print("OC{$ORD,$DIM}: 2/3")
    Base.show(io::IO,::OperatorConstructor{ORD,DIM,_1,_2,_3,_4,TBK,Union{},TB2}) where {ORD,DIM,_1,_2,_3,_4,TBK,TB2} = print("OC{$ORD,$DIM}: 2/3")
    Base.show(io::IO,::OperatorConstructor{ORD,DIM,_1,_2,_3,_4,Union{},TB1,TB2}) where {ORD,DIM,_1,_2,_3,_4,TB1,TB2} = print("OC{$ORD,$DIM}: 2/3")
    Base.show(io::IO,::OperatorConstructor{ORD,DIM}) where {ORD,DIM} = print("OC{$ORD,$DIM}: 3/3")
end

function get_vals(OC::OperatorConstructor{2,DIM,CS,TBL1,TBL2},args::Vararg{Number,M}) where {M,DIM,CS,TD,TBL1,TBL2}
    N = OC.N

    h = 1 .+ OC.definition.hd[DIM]/args[1]
    vals = (1 ./h)/OC.dx[DIM]^2
    bulk_on_site = -(vals[2:end-2] + vals[3:end-1])
    bnd1_on_site = -(vals[1] + vals[2])
    bnd2_on_site = -(vals[end-1] + vals[end])
    bulk_nn = vals[2:N[DIM]]
    return vcat(bulk_on_site,bulk_nn,bulk_nn,bnd1_on_site,bnd2_on_site)
end
function get_vals(OC::OperatorConstructor{2,1,Polar,TBL1,TBL2},args::Vararg{Number,M}) where {M,TD,TBL1,TBL2}
    DIM=1

    N = OC.N

    r = OC.definition.xd[DIM]
    rf = r .+ OC.definition.fd[DIM]/args[1]
    h = 1 .+ OC.definition.hd[DIM]/args[1]

    vals = (rf./h)/OC.dx[DIM]^2
    bulk_on_site = -(vals[2:end-2] + vals[3:end-1])
    bnd1_on_site = -(vals[1] + vals[2])
    bnd2_on_site = -(vals[end-1] + vals[end])
    bulk_nn = vals[2:N[DIM]]
    return vcat(bulk_on_site,bulk_nn,bulk_nn,bnd1_on_site,bnd2_on_site)
end


function get_vals(OC::OperatorConstructor{2,DIM,CS,TBL1,TBL2},BC::T,args...) where {DIM,CS,T<:AbstractLocalBC{DIM,SIDE},TD,TBL1,TBL2} where SIDE
    @assert (TBL1<:noBL && TBL2<:noBL) || length(args)>0 "pml requires frequency argument. Try with additional argument"
    weights = BC.weights
    if CS<:Polar && DIM==1
        fd = SIDE==1 ? OC.definition.fd[DIM][1] : OC.definition.fd[DIM][end]
        rf = OC.Δ[DIM][SIDE] .+ fd/args[1]
    else
        rf = 1
    end
    ind = SIDE==1 ? 1 : OC.N[DIM]
    hd = 1 .+ OC.definition.hd[DIM][ind]/args[1]
    return [(rf./hd).*weights[1]/OC.dx[DIM]^2]
end

function get_vals(OC::OperatorConstructor{2,DIM,CS},BC::PeriodicBC{DIM,SIDE},k::Number,ka::Number,kb::Number,ind::Int) where {DIM,SIDE,CS}
    @assert !(CS<:Polar && DIM==1) "combination DIM=1, CS=Polar, BC=Periodic not recognized"
    perp_ind = mod1(ind+1,2)
    I1, I2, J1, J2 = Int[], Int[], Int[], Int[]
    vals = ComplexF64[]
    ns = Int[]
    for i ∈ eachindex(BC.shifts[perp_ind])
        I1 = vcat(I1,BC.I1[perp_ind][i])
        I2 = vcat(I2,BC.I2[perp_ind][i])
        J1 = vcat(J1,BC.J1[perp_ind][i])
        J2 = vcat(J2,BC.J2[perp_ind][i])
        vals = vcat(vals,BC.weights[perp_ind][i])
        ns = vcat(ns,fill(BC.shifts[perp_ind][i],length(BC.weights[perp_ind][i])))
    end
    (K,κ) = ind==1 ? (ka,kb) : (kb,ka)
    (R,r) = ind==1 ? (BC.lattice.a,BC.lattice.b) : (BC.lattice.b,BC.lattice.a)
    if !isinf(r)
        ϕ=κ*r
    else
        ϕ=0.0
    end
    idx = 1
    weights = Array{Array{ComplexF64,1},1}(undef,length(BC.shifts[ind]))
    for i ∈ eachindex(BC.shifts[ind])
        rng = idx:(idx+length(BC.weights[ind][i])-1)
        vals[rng] .*= exp.(-1im*ns[rng]*ϕ)
        weights[i] = vals[rng]
        idx += length(BC.weights[ind][i])
    end
    h = SIDE==1 ? 1 .+ OC.definition.hd[DIM][1]/k : 1 .+ OC.definition.hd[DIM][end]/k
    for i ∈ eachindex(weights)
        weights[i] .*= (1/h)/OC.dx[DIM]^2
    end
    return weights
end

function get_vals(OC::OperatorConstructor{2,DIM,CS},BC::MatchedBC{DIM,SIDE},args...) where {DIM,CS,SIDE}
    vals = Array{Array{ComplexF64,1},1}(undef,length(BC.QNs))
    weights = BC.weights
    R = (CS<:Polar && DIM==1) ? OC.Δ[DIM][SIDE] : 1
    for i ∈ eachindex(vals)
        vals[i] = R.*weights[i]/OC.dx[DIM]^2
    end
    return vals
end

#######################################################################################

struct DifferentialOperator{TD,TS}
    N::NTuple{2,Int}
    d::Array{RowColVal{TD},1}
    f::Array{Function,1}
    s::RowColVal{TS}

    DifferentialOperator(d::RowColVal{TD},s::RowColVal{TS},f::TF) where {TD,TS,TF<:Function} = new{TD,TS}(d.N,[d],[f],s)
    DifferentialOperator(d::Array{RowColVal{TD},1},f::Array{Function,1},s::RowColVal{TS}) where {TD,TS} = new{TD,TS}(d[1].N,d,f,s)

    Base.show(io::IO,DO::DifferentialOperator) = print("Differential Operator")

    Base.vcat(D::DifferentialOperator) = D
    function Base.vcat(D1::DifferentialOperator{TD1,TS1},D2::DifferentialOperator{TD2,TS2}) where {TD1,TD2,TS1,TS2}
        @assert D1.N==D2.N "cannot combine Differential Operators with different N, here D1.N=$(D1.N), D2.N=$(D2.N)"
        d = vcat(D1.d,D2.d)
        f = vcat(D1.f,D2.f)
        s = D1.s
        new{promote_type(TD1,TD2),promote_type(TS1,TS2)}(D1.N,d,f,s)
    end
    Base.vcat(D::DifferentialOperator,Ds...) = vcat(D,vcat(Ds...))

    SparseArrays.sparse(DO::DifferentialOperator) = (sparse.(DO.d),DO.f,sparse(DO.s))
    function SparseArrays.sparse(D1::DifferentialOperator,D2::DifferentialOperator)
        d1,_,s1 = sparse(D1)
        d2,_,s2 = sparse(D2)
        T = promote_type(eltype(d1[1]),eltype(d2[1]))
        D = Array{SparseMatrixCSC{T,Int64},1}(undef,length(d1)+length(d2))
        for i ∈ eachindex(d1)
            D[i] = s2*d1[i]
        end
        for i ∈ eachindex(d2)
            D[length(d1) + i] = s1*d2[i]
        end
        return (D,vcat(D1.f,D2.f),s1*s2)
    end
end


"""
    diff_op(definition, args...)
"""
function diff_op(OD::OperatorDefinition, args...)
    if !(typeof(OD.bls)<:NTuple{2,NTuple{2,noBL}})
        @assert length(args)≥1 "must provide at least one argument for PMLs: k"
    elseif length(args)==0
        args = (1,)
    else
        nothing
    end
    OC = OperatorConstructor(OD)
    blk_op = bulk_op(OD,OC,args...)
    bnd_op = boundary_op(OD,OC,args...)
    return vcat(blk_op,bnd_op)
end


"""
    bulk_op(definition, args...)
"""
bulk_op(OD::OperatorDefinition, args...) = bulk_op(OD, OperatorConstructor(OD), args...)
function bulk_op(OD::OperatorDefinition{ORD,DIM,CS}, OC::OperatorConstructor, args...) where {ORD,DIM,CS}
    OC = OC(args...) # populate bulk nodes

    h = 1 .+ OC.definition.h[DIM]/args[1]
    if CS<:Polar && DIM==1
        rf = OC.definition.x[DIM] + OC.definition.f[DIM]/args[1]
        s = RowColVal(OD.N,rf.*h,DIM)
    else
        s = RowColVal(OD.N,h,DIM)
    end

    return DifferentialOperator(OC.bulk,s,(a...)->1)
end


"""
    boundary_op(definition, args...)
"""
boundary_op(OD::OperatorDefinition, args...) = boundary_op(OD,OperatorConstructor(OD),args...)
function boundary_op(OD::OperatorDefinition{ORD,DIM,CS},OC::OperatorConstructor, args...) where {ORD,DIM,CS}
    OC = OC(OD.bcs[DIM][1],args...)
    OC = OC(OD.bcs[DIM][2],args...)
    rcv = vcat(OC.bnd1,OC.bnd2)

    if typeof(OD.bcs[DIM][1])<:PeriodicBC
        if args[4]==0
            fs1 = map(x->((k,_...)->OD.bcs[DIM][1](k,args[1],args[2],1)[x]),1:length(OC.bnd1))
            fs2 = map(x->((k,_...)->OD.bcs[DIM][2](k,args[1],args[2],1)[x]),1:length(OC.bnd2))
        elseif args[4]==1
            fs1 = map(x->((k,κ)->OD.bcs[DIM][1](k,κ,args[2],1)[x]),1:length(OC.bnd1))
            fs2 = map(x->((k,κ)->OD.bcs[DIM][2](k,κ,args[2],1)[x]),1:length(OC.bnd2))
        elseif args[4]==2
            fs1 = map(x->((k,κ)->OD.bcs[DIM][1](k,args[3],κ,2)[x]),1:length(OC.bnd1))
            fs2 = map(x->((k,κ)->OD.bcs[DIM][2](k,args[3],κ,2)[x]),1:length(OC.bnd2))
        end
    else
        fs1 = map(x->((k,_...)->OD.bcs[DIM][1](k)[x]),1:length(OC.bnd1))
        fs2 = map(x->((k,_...)->OD.bcs[DIM][2](k)[x]),1:length(OC.bnd2))
    end
    fs = vcat(fs1,fs2)

    h = 1 .+ OC.definition.h[DIM]/args[1]
    if CS<:Polar && DIM==1
        rf = OC.definition.x[DIM] + OC.definition.f[DIM]/args[1]
        s = RowColVal(OD.N,rf.*h,DIM)
    else
        s = RowColVal(OD.N,h,DIM)
    end

    return DifferentialOperator(rcv, fs, s)
end


end # module
