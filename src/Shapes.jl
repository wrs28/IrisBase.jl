module Shapes

using Formatting,
RecipesBase

export AbstractShape,
Circle,
Ellipse,
Square,
Rectangle,
Parallelogram,
Universe

abstract type AbstractShape end


"""
    Circle(R, x0, y0)
"""
struct Circle <: AbstractShape
    R::Float64
    x0::Float64
    y0::Float64

    Circle(R::Number,x0::Number,y0::Number) = new(R,x0,y0)
    (c::Circle)(x,y) = hypot((x-c.x0), (y-c.y0)) ≤ c.R
end
Base.show(io::IO, circle::Circle) = begin
    print(io, "Circle(R=$(fmt("2.2f",circle.R)), x0=$(fmt("2.2f",circle.x0)), y0=$(fmt("2.2f",circle.y0)))")
end
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


"""
    Ellipse(a, b, x0, y0, θ)
"""
struct Ellipse <: AbstractShape
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
end
Base.show(io::IO, ellipse::Ellipse) = begin
    print(io, "Ellipse(a=$(fmt("2.2f",ellipse.a)), b=$(fmt("2.2f",ellipse.b)), x0=$(fmt("2.2f",ellipse.x0)), y0=$(fmt("2.2f",ellipse.y0)), θ=$(fmt("2.2f",ellipse.θ)))")
end
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


"""
    Square(a,x0,y0,θ)
"""
struct Square <: AbstractShape
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
end
Base.show(io::IO, square::Square) = begin
    print(io, "Square(a=$(fmt("2.2f",square.a)), x0=$(fmt("2.2f",square.x0)), y0=$(fmt("2.2f",square.y0)), θ=$(fmt("2.2f",square.θ)))")
end
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


"""
    Rectangle(a, b, x0, y0, θ)
"""
struct Rectangle <: AbstractShape
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
end
Base.show(io::IO, rect::Rectangle) = begin
    print(io, "Rectangle(a=$(fmt("2.2f",rect.a)), b=$(fmt("2.2f",rect.b)), x0=$(fmt("2.2f",rect.x0)), y0=$(fmt("2.2f",rect.y0)), θ=$(fmt("2.2f",rect.θ)))")
end
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


"""
    Parallelogram(a, b, α, x0, y0, θ)
"""
struct Parallelogram <: AbstractShape
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
end
Base.show(io::IO, par::Parallelogram) = begin
    print(io, "Parallelogram(a=$(fmt("2.2f",par.a)), b=$(fmt("2.2f",par.b)), α=$(fmt("2.2f",par.α)), x0=$(fmt("2.2f",par.x0)), y0=$(fmt("2.2f",par.y0)), θ=$(fmt("2.2f",par.θ)))")
end
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


struct Universe <: AbstractShape
    Universe() = new()
    (::Universe)(x,y) = true
end
Base.show(io::IO, ::Universe) = print(io, "Universe()")


"""
    rotate(x,y,cosθ,sinθ)
"""
function rotate(x::Number,y::Number,cosθ::Number,sinθ::Number)
    if !iszero(sinθ)
        xrot = +cosθ*x - sinθ*y
        yrot = +sinθ*x + cosθ*y
    else
        xrot = x
        yrot = y
    end
    return xrot, yrot
end
function rotate(x,y,cosθ::Number,sinθ::Number)
    xyrot = rotate.(x,y,cosθ,sinθ)
    xrot = Array{Float64}(undef,size(xyrot))
    yrot = Array{Float64}(undef,size(xyrot))
    for i ∈ eachindex(xyrot)
        xrot[i] = xyrot[i][1]
        yrot[i] = xyrot[i][2]
    end
    return xrot, yrot
end


end #module
