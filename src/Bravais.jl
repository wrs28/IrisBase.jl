module Bravais

using Formatting,
RecipesBase

export BravaisLattice,
bravais_coordinates_unit_cell,
bravais_coordinates_unit_cell!,
bravais_coordinates

"""
    lattice = BravaisLattice(; a=Inf, α=0, b=Inf, β=π/2, x0=0, y0=0)

bravais lattice with:

- `a` is the length of the first vector (by default aligned along x-axis).

- `b` is length of second vector.

- `α` is the angle of first primitive vector.

- `β` is the angle of second primitive vector

- `x0`, `y0` define the origin of the lattice.
"""
struct BravaisLattice
    a::Float64
    b::Float64
    α::Float64
    β::Float64
    x0::Float64
    y0::Float64

    sin²θ::Float64
    cosθ::Float64
    R::Array{Float64,2}#SArray{Tuple{2,2},Float64,2,4}#
    Rᵀ::Array{Float64,2}#SArray{Tuple{2,2},Float64,2,4}#
    v1::Array{Float64,1}#SArray{Tuple{2},Float64,1,2}#
    v2::Array{Float64,1}#SArray{Tuple{2},Float64,1,2}#

    function BravaisLattice(;a::Number=Inf, α::Number=0, b::Number=Inf, β::Number=π/2, x0::Number=0, y0::Number=0)
        θ = β-α
        if iszero(θ)
            throw(ArgumentError("lattice angle θ=$(θ) cannot be zero"))
        end
        sinθ, cosθ = sincos(θ)
        # R = @SMatrix [ cos(α) -sin(α); sin(α)  cos(α)]
        R = [ cos(α) -sin(α); sin(α)  cos(α)]
        Rᵀ = transpose(R)
        # v1 = R*SVector(1,0)
        v1 = R*[1,0]
        # v2 = R*SVector(cosθ,sinθ)
        v2 = R*[cosθ,sinθ]
        new(float(a), float(b), float(α), float(β), float(x0), float(y0),
            sinθ^2, cosθ, R, Rᵀ, v1, v2)
    end

    function Base.show(io::IO, bvl::BravaisLattice)
        if !get(io, :sub, false)
            print(io, typeof(bvl), ": \n",
            "\tprimitive vector 1: ", fmt("3.2f",bvl.a), ", ∠", fmt("3.2f",(mod2pi(bvl.α))*180/pi), "°\n",
            "\tprimitive vector 2: ", fmt("3.2f",bvl.b), ", ∠", fmt("3.2f",(mod2pi(bvl.β))*180/pi), "°\n",
            "\torigin: (", fmt("3.2f",bvl.x0), ", ", fmt("3.2f",bvl.y0), ")")
        elseif !get(io, :sub1, false)
            print(io,
            "\n\tprimitive vector 1: ", fmt("3.2f",bvl.a), ", ∠", fmt("3.2f",(mod2pi(bvl.α))*180/pi), "°\n",
            "\tprimitive vector 2: ", fmt("3.2f",bvl.b), ", ∠", fmt("3.2f",(mod2pi(bvl.β))*180/pi), "°\n",
            "\torigin: (", fmt("3.2f",bvl.x0), ", ", fmt("3.2f",bvl.y0), ")")
        else
            print(io,
            "\n\t\tprimitive vector 1: ", fmt("3.2f",bvl.a), ", ∠", fmt("3.2f",(mod2pi(bvl.α))*180/pi), "°\n",
            "\t\tprimitive vector 2: ", fmt("3.2f",bvl.b), ", ∠", fmt("3.2f",(mod2pi(bvl.β))*180/pi), "°\n",
            "\t\torigin: (", fmt("3.2f",bvl.x0), ", ", fmt("3.2f",bvl.y0), ")")
        end
    end
end # end of struct BravaisLatice

"""
    lattice = BravaisLattice(lattice; :key1 => value1, :key2 => value2, ...)

new lattice from old lattice with modified fields.
"""
BravaisLattice(lattice::BravaisLattice; a=lattice.a, b=lattice.b, α=lattice.α, β=lattice.β, x0=lattice.x0, y0=lattice.y0) = BravaisLattice(;a=a, b=b, α=α, β=β, x0=x0, y0=y0)

################################################################################
# COORDINATE TRANSFORMATIONS
################################################################################
"""
    xb, yb = bravais_coordinates_unit_cell(x, y, lattice)

maps cartesian (`x`,`y`) into cartesian (`xb`,`yb`) unit cell specified by
lattice
"""
function bravais_coordinates_unit_cell(x::Real, y::Real, lattice::BravaisLattice)
    a, b, α, β, v1, v2, x0, y0 = lattice.a, lattice.b, lattice.α, lattice.β, lattice.v1, lattice.v2, lattice.x0, lattice.y0
    if (isinf(a) && iszero(β-π/2)) || (isinf(b) && iszero(α)) || (isinf(a) && isinf(b))
        return float(x), float(y)
    end

    p1, p2 = bravais_coordinates(x, y, lattice)
    p1 = mod(p1, a)
    p2 = mod(p2, b)

    xb = v1[1]*p1 + v2[1]*p2 + x0
    yb = v1[2]*p1 + v2[2]*p2 + y0

    return xb, yb
end

function bravais_coordinates_unit_cell(x::AbstractArray, y::AbstractArray, lattice::Union{Tlat,AbstractArray{Tlat}}) where Tlat<:BravaisLattice
    xb = Array{Float64}(undef, size(x)...)
    yb = Array{Float64}(undef, size(y)...)
    bravais_coordinates_unit_cell!(xb,yb,x,y,lattice)
    return xb, yb
end

function bravais_coordinates_unit_cell!(xb::AbstractArray, yb::AbstractArray, x::AbstractArray, y::AbstractArray, lattice::BravaisLattice)
    for k ∈ eachindex(x)
        xb[k], yb[k] = bravais_coordinates_unit_cell(x[k],y[k],lattice)
    end
    return nothing
end
function bravais_coordinates_unit_cell!(xb::AbstractArray, yb::AbstractArray, x::AbstractArray, y::AbstractArray, lattices::AbstractArray{Tlat}) where Tlat<:BravaisLattice
    for k ∈ eachindex(x)
        xb[k], yb[k] = bravais_coordinates_unit_cell(x[k],y[k],lattices[k])
    end
    return nothing
end

"""
    p1, p2 = bravais_coordinates(x,y,lattice)

coordinates in bravais frame (i.e. (x,y) = p1*v1 + p2*v2)
"""
function bravais_coordinates(x::Real, y::Real, lattice::BravaisLattice)
    a, b, α, β, v1, v2, x0, y0 = lattice.a, lattice.b, lattice.α, lattice.β, lattice.v1, lattice.v2, lattice.x0, lattice.y0
    if (isinf(a) && iszero(β-π/2)) || (isinf(b) && iszero(α)) || (isinf(a) && isinf(b))
        return float(x-x0), float(y-y0)
    end

    x += -x0
    y += -y0

    rv1 = v1[1]*x + v1[2]*y
    rv2 = v2[1]*x + v2[2]*y

    p1 = (rv1 - rv2*lattice.cosθ)/lattice.sin²θ
    p2 = (rv2 - rv1*lattice.cosθ)/lattice.sin²θ
    return p1, p2
end
function bravais_coordinates(x::AbstractArray,y::AbstractArray,lattice::BravaisLattice)
    p1 = Array{Float64,2}(undef,length(x),length(y))
    p2 = Array{Float64,2}(undef,length(x),length(y))
    bravais_coordinates!(p1,p2,x,y,lattice)
    return p1, p2
end
function bravais_coordinates!(p1::AbstractArray,p2::AbstractArray,x::AbstractArray,y::AbstractArray,lattice::BravaisLattice)
    if size(x)!==size(y)
        for i ∈ eachindex(x), j ∈ eachindex(y)
            p1[i,j],p2[i,j] = bravais_coordinates(x[i],y[j],lattice)
        end
    else
        for i ∈ eachindex(x)
            p1[i],p2[i] = bravais_coordinates(x[i],y[i],lattice)
        end
    end
    return nothing
end
bravais_coordinates(x::Real, y::AbstractArray, args...) = bravais_coordinates(hcat(x),y,args...)
bravais_coordinates(x::AbstractArray, y::Real, args...) = bravais_coordinates(x,vcat(y),args...)

########################################################################
### PLOTTING
########################################################################
"""
    plot(lattice::BravaisLattice, N=[4,4])

scatter plot of lattice vectors from -`N[i]`:`+N[i]` for each dimension i.

If `N` is a scalar, uses same span for both dimensions.
"""
@recipe function f(lattice::BravaisLattice, N::Array{Int,1}=[4,4])
    seriestype --> :scatter
    legend --> false
    aspect_ratio --> 1
    xlims --> [-N[1]*(lattice.a*lattice.v1[1]+lattice.b*lattice.v2[1]), +N[1]*(lattice.a*lattice.v1[1]+lattice.b*lattice.v2[1])]/2
    ylims --> [-N[2]*(lattice.a*lattice.v1[2]+lattice.b*lattice.v2[2]), +N[2]*(lattice.a*lattice.v1[2]+lattice.b*lattice.v2[2])]/2
    @series begin
        x, y = Float64[], Float64[]
        for n ∈ -N[1]:N[1], m ∈ -N[2]:N[2]
            temp = n*lattice.a*lattice.v1 + m*lattice.b*lattice.v2
            push!(x,temp[1]); push!(y,temp[2])
        end
        (x, y)
    end
    @series begin
        x, y = Float64[], Float64[]
        for i ∈ 1:5
            n = [0,1,1,0,0][i]
            m = [0,0,1,1,0][i]
            temp = n*lattice.a*lattice.v1 + m*lattice.b*lattice.v2
            push!(x,temp[1]); push!(y,temp[2])
        end
        seriestype --> :path
        fill --> (0,.15, :red)
        (x, y)
    end
end
@recipe f(lattice::BravaisLattice, N::Int) = (lattice, [N,N])

end # module
