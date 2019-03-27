module DielectricFunctions

export piecewise_constant_ε,
piecewise_constant_F

"""
    piecewise_constant_ε(x,y,params)
"""
function piecewise_constant_ε(x,y,params::Dict)
    n2 = get!(params,:n₂,0)
    n1 = get!(params,:n₁,1)
    return complex(n1,n2)^2
end


"""
    piecewise_constant_F(x,y,params)
"""
piecewise_constant_F(x,y,params) = get!(params,:F,0.0)

end # module
