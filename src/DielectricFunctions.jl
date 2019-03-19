module DielectricFunctions

export piecewise_constant_ε,
piecewise_constant_F


"""
    piecewise_constant_ε(x,y,params)
"""
function piecewise_constant_ε(x,y,params::Dict,n1::Array{Float64,1},n2::Array{Float64,1})

    n2[1] = get!(params,:n₂,0)
    n2[1] = get!(params,:n2,n2[1])

    n1[1] = get!(params,:n₁,1)
    n1[1] = get!(params,:n1,n1[1])

    return complex(n1[1],n2[1])^2
end


"""
    piecewise_constant_F(x,y,params)
"""
function piecewise_constant_F(x,y,params::Dict,F::Array{Float64,1})
    F[1] = get!(params,:F,0.0)
    return F[1]
end


end # module
