module Shapes

using Formatting,
RecipesBase

export AbstractShape,
Circle,
Ellipse,
Square,
Rectangle,
Parallelogram,
DeformedDisk,
Universe

abstract type AbstractShape end
abstract type AbstractParallelogram <: AbstractShape end
abstract type AbstractRectangle <: AbstractParallelogram end
abstract type AbstractSquare <: AbstractRectangle end
abstract type AbstractEllipse <: AbstractShape end
abstract type AbstractCircle <: AbstractEllipse end

"""
    Circle(R, x0, y0)
"""
struct Circle <: AbstractCircle
    R::Float64
    x0::Float64
    y0::Float64

    Circle(R::Number,x0::Number,y0::Number) = new(R,x0,y0)
    (c::Circle)(x,y) = hypot((x-c.x0), (y-c.y0)) ≤ c.R

    Base.show(io::IO, circle::Circle) = print(io, "Circle(R=$(fmt("2.2f",circle.R)), x0=$(fmt("2.2f",circle.x0)), y0=$(fmt("2.2f",circle.y0)))")

    @recipe function f(cr::Circle)
        θ = LinRange(0,2π,101)
        x = @. cr.x0 + cr.R*cos(θ)
        y = @. cr.x0 + cr.R*sin(θ)
        seriestype --> :path
        fill --> (0,.15,:red)
        aspect_ratio --> 1
        legend-->false
        (x,y)
    end
end

"""
    Ellipse(a, b, x0, y0, θ)
"""
struct Ellipse <: AbstractEllipse
    a::Float64
    b::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    cosθ::Float64
    sinθ::Float64

    Ellipse(a::Number,b::Number,x0::Number,y0::Number,θ::Number) = new(float(a),float(b),float(x0),float(y0),float(θ),cos(θ),sin(θ))
    function (e::Ellipse)(x,y)
        xrot, yrot = rotate(x-e.x0, y-e.y0, e.cosθ, e.sinθ)
        return hypot(xrot/e.a, yrot/e.b) ≤ 1
    end

    Base.show(io::IO, ellipse::Ellipse) = print(io, "Ellipse(a=$(fmt("2.2f",ellipse.a)), b=$(fmt("2.2f",ellipse.b)), x0=$(fmt("2.2f",ellipse.x0)), y0=$(fmt("2.2f",ellipse.y0)), θ=$(fmt("2.2f",ellipse.θ)))")

    @recipe function f(el::Ellipse)
        θ = LinRange(0,2π,101)
        x = @. el.a*cos(θ)
        y = @. el.b*sin(θ)
        seriestype --> :path
        fill --> (0,.15,:red)
        aspect_ratio --> 1
        legend-->false
        xrot, yrot = rotate(x,y,el.cosθ,el.sinθ)
        el.x0 .+ xrot, el.y0 .+ yrot
    end
end


"""
    Square(a,x0,y0,θ)
"""
struct Square <: AbstractSquare
    a::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    cosθ::Float64
    sinθ::Float64

    function Square(a::Number,x0::Number,y0::Number,θ::Number)
        new(float(a),float(x0),float(y0),float(θ),cos(θ),sin(θ))
    end
    function (s::Square)(x,y)
        xrot, yrot = rotate(x-s.x0, y-s.y0, s.cosθ, s.sinθ)
        return 0 ≤ xrot ≤ s.a  &&  0 ≤ yrot ≤ s.a
    end

    Base.show(io::IO, square::Square) = print(io, "Square(a=$(fmt("2.2f",square.a)), x0=$(fmt("2.2f",square.x0)), y0=$(fmt("2.2f",square.y0)), θ=$(fmt("2.2f",square.θ)))")

    @recipe function f(sq::Square)
        x = sq.a*cumsum([0, +1,  0, -1,  0])
        y = sq.a*cumsum([0,  0, +1,  0, -1])
        seriestype --> :path
        fill --> (0,.15,:red)
        aspect_ratio --> 1
        legend --> false
        xrot, yrot = rotate(x,y,sq.cosθ,sq.sinθ)
        sq.x0 .+ xrot, sq.y0 .+ yrot
    end
end


"""
    Rectangle(a, b, x0, y0, θ)
"""
struct Rectangle <: AbstractRectangle
    a::Float64
    b::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    cosθ::Float64
    sinθ::Float64

    function Rectangle(a::Number,b::Number,x0::Number,y0::Number,θ::Number)
        new(float(a),float(b),float(x0),float(y0),float(θ),cos(θ),sin(θ))
    end

    function (r::Rectangle)(x,y)
        xrot, yrot = rotate(x-r.x0, y-r.y0, r.cosθ, r.sinθ)
        return 0 ≤ xrot ≤ r.a  &&  0 ≤ yrot ≤ r.a
    end

    Base.show(io::IO, rect::Rectangle) = print(io, "Rectangle(a=$(fmt("2.2f",rect.a)), b=$(fmt("2.2f",rect.b)), x0=$(fmt("2.2f",rect.x0)), y0=$(fmt("2.2f",rect.y0)), θ=$(fmt("2.2f",rect.θ)))")

    @recipe function f(rc::Rectangle)
        x = cumsum([0, rc.a,      0, -rc.a, 0])
        y = cumsum([0,     0, +rc.b,     0, -rc.b])
        seriestype --> :path
        fill --> (0,.15,:red)
        aspect_ratio --> 1
        legend --> false
        xrot, yrot = rotate(x,y,rc.cosθ,rc.sinθ)
        rc.x0 .+ xrot, rc.y0 .+ yrot
    end
end


"""
    Parallelogram(a, b, α, x0, y0, θ)
"""
struct Parallelogram <: AbstractParallelogram
    a::Float64
    b::Float64
    α::Float64
    x0::Float64
    y0::Float64
    θ::Float64
    cosθ::Float64
    sinθ::Float64
    tanα::Float64
    cosα::Float64

    function Parallelogram(a::Number, b::Number, α::Number, x0::Number, y0::Number, θ::Number)
        new(float(a),float(b),float(α),float(x0),float(y0),float(θ),cos(θ),sin(θ),tan(α),cos(α))
    end

    function (p::Parallelogram)(x,y)
        xrot, yrot = rotate(x-p.x0, y-p.y0, p.cosθ, p.sinθ)
        return p.tanα*(xrot-p.a) ≤ yrot ≤ p.tanα*xrot  &&  0 ≤ yrot ≤ p.cosα*p.b
    end

    Base.show(io::IO, par::Parallelogram) = print(io, "Parallelogram(a=$(fmt("2.2f",par.a)), b=$(fmt("2.2f",par.b)), α=$(fmt("2.2f",par.α)), x0=$(fmt("2.2f",par.x0)), y0=$(fmt("2.2f",par.y0)), θ=$(fmt("2.2f",par.θ)))")

    @recipe function f(pg::Parallelogram)
        x = cumsum([0, pg.a, +pg.b*cos(pg.α), -pg.a, -pg.b*cos(pg.α)])
        y = cumsum([0, 0, +pg.b*sin(pg.α), 0, -pg.b*sin(pg.α)])
        seriestype --> :path
        fill --> (0,.15,:red)
        aspect_ratio --> 1
        legend --> false
        xrot, yrot = rotate(x,y,pg.cosθ,pg.sinθ)
        pg.x0 .+ xrot, pg.y0 .+ yrot
    end
end


"""
    DeformedDisk(R, x0, y0, M, a, φ)
`R` is radius, `x0` and `y0` are origin, `M`
is array of multipole integers, `a` is array of
amplitudes, `φ` is array of angles
"""
struct DeformedDisk{N} <: AbstractShape
    R::Float64
    M::Vector{Int}
    a::Vector{Float64}
    φ::Vector{Float64}
    x0::Float64
    y0::Float64

    function DeformedDisk(R::Number,x0::Number,y0::Number,M::AbstractArray,a,φ)
        new{length(M)}(float(R),M,a,φ,float(x0),float(y0))
    end
    function (d::DeformedDisk)(x,y)
        θ = atan(y-d.y0,x-d.x0)
        return hypot(x,y) ≤ d.R + sum(d.a.*(map((x,y)->cos(x*(θ-y)),d.M,d.φ)))
    end

    Base.show(io::IO, d::DeformedDisk) = print(io,"DeformedDisk")#(a=$(fmt("2.2f",ellipse.a)), b=$(fmt("2.2f",ellipse.b)), x0=$(fmt("2.2f",ellipse.x0)), y0=$(fmt("2.2f",ellipse.y0)), θ=$(fmt("2.2f",ellipse.θ)))")

    @recipe function f(d::DeformedDisk)
        θ = LinRange(0,2π,101)
        r = d.R .+ sum(d.a.*(map((x,y)->cos.(x*(θ.-y)),d.M,d.φ)))
        x = @. r*cos(θ)
        y = @. r*sin(θ)
        seriestype --> :path
        fill --> (0,.15,:red)
        aspect_ratio --> 1
        legend-->false
        d.x0 .+ x, d.y0 .+ y
    end
end


"""
    Universe()
"""
struct Universe <: AbstractShape
    Universe() = new()
    (::Universe)(x,y) = true

    Base.show(io::IO, ::Universe) = print(io, "Universe()")
end


"""
    rotate(x,y,cosθ,sinθ)
"""
function rotate(x::Real,y::Real,cosθ::Real,sinθ::Real)
    (xrot,yrot) = iszero(sinθ) ? (x,y) : (cosθ*x-sinθ*y, sinθ*x+cosθ*y)
    return xrot, yrot
end

# function rotate(x,y,cosθ::Number,sinθ::Number)
#     xyrot = rotate.(x,y,cosθ,sinθ)
#     xrot = Array{Float64}(undef,size(xyrot))
#     yrot = Array{Float64}(undef,size(xyrot))
#     for i ∈ eachindex(xyrot)
#         xrot[i] = xyrot[i][1]
#         yrot[i] = xyrot[i][2]
#     end
#     return xrot, yrot
# end

end #module
