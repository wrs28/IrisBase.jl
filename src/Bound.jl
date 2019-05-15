module Bound

using ..Shapes

struct Bounds{TS,TB}
    shape::TS
    bcs::TB
end # struct


function unit_normal!(Nx,Ny,C::AbstractCircle,x,y)
    for i ∈ eachindex(Nx)
        Nx[i], Ny[i] = unit_normal(C,x[i],y[i])
    end
    return nothing
end
function unit_normal(C::AbstractCircle,x,y)
    r = hypot(x-C.x0,y-C.y0)
    return x/r,y/r
end

end # module
