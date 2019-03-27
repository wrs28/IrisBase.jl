# TODO: sub_pixel_smoothing in 2-dim case, take care of smoothing along edges
# TODO: fix Bravais in 1-dimension

################################################################################
###############   SUB PIXEL SMOOTHING    #######################################
################################################################################
function construct_εFr(bnd::Boundary, dis::Discretization, sys::System)
    x = dis.X
    y = dis.Y
    domains = Array{Int}(undef,size(x)...)
    which_domains!(domains, x, y, bnd, sys)
    r = domains

    xb, yb = Array{Float64}(undef,size(x)...), Array{Float64}(undef,size(y)...)
    bravais_coordinates_unit_cell!(xb,yb,x,y,domains,sys)

    # n1, n2, f = 1.0, 0.0, 0.0
    ε = Array{ComplexF64}(undef,dis.N[1],dis.N[2])
    F = Array{Float64}(undef,dis.N[1],dis.N[2])
    for i ∈ eachindex(r)
        j = r[i]
        ε[i] = sys.ε_by_region[j](x[i],y[i],sys.params_by_region[j])
        F[i] = sys.F_by_region[j](x[i],y[i],sys.params_by_region[j])
    end
    return ε, F, r
end


"""
    ε, regions = sub_pixel_smoothing(bnd, dis, sys; display=false)
"""
function sub_pixel_smoothing!(bnd::Boundary, dis::Discretization{Cartesian}, sys::System, ε, F, r; display::Bool=false)
    x = dis.X
    y = dis.Y
    sub_pixel_num = dis.sub_pixel_num
    if dis.N[2] == 1

    elseif dis.N[1] == 1

    else
        sub_x = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        sub_y = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        xb = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        yb = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        domains = Array{Int}(undef, sub_pixel_num, sub_pixel_num)
        ε_temp = Array{ComplexF64}(undef,sub_pixel_num, sub_pixel_num)
        F_temp = Array{Float64}(undef,sub_pixel_num, sub_pixel_num)
        if display
            pg = Progress((size(x,1)-2)*(size(y,2)-2), PROGRESS_UPDATE_TIME::Float64, "sub-pixel smoothing ")
        end
        for h ∈ CartesianIndices(x)
            i,j = h[1],h[2]
            if all((1,1).<(i,j).<size(x))
                nearestNeighborFlag = r[i,j]!==r[i,j+1] || r[i,j]!==r[i,j-1] || r[i,j]!==r[i+1,j] || r[i,j]!==r[i-1,j]
                nextNearestNeighborFlag = r[i,j]!==r[i+1,j+1] || r[i,j]!==r[i-1,j-1] || r[i,j]!==r[i+1,j-1] || r[i,j]!==r[i-1,j+1]
                if nearestNeighborFlag || nextNearestNeighborFlag
                    x_min = (x[i,j]+x[i-1,j])/2
                    x_max = (x[i,j]+x[i+1,j])/2
                    y_min = (y[i,j]+y[i,j-1])/2
                    y_max = (y[i,j]+y[i,j+1])/2
                    for k ∈ CartesianIndices(sub_x)
                        sub_x[k] = x_min + (k[1]-1)*(x_max-x_min)/(sub_pixel_num-1)
                        sub_y[k] = y_min + (k[2]-1)*(y_max-y_min)/(sub_pixel_num-1)
                    end
                    which_domains!(domains, sub_x, sub_y, bnd, sys)
                    bravais_coordinates_unit_cell!(xb, yb, sub_x, sub_y, domains, sys)
                    for l ∈ eachindex(xb)
                        d = domains[l]
                        ε_temp[l] = sys.ε_by_region[d](sub_x[d],sub_y[d],sys.params_by_region[d])
                        F_temp[l] = sys.F_by_region[d](sub_x[d],sub_y[d],sys.params_by_region[d])
                    end
                    ε[i,j] = mean(ε_temp)
                    F[i,j] = mean(F_temp)
                end
                if display
                    next!(pg)
                end
            end
        end
    end
    return nothing
end
function sub_pixel_smoothing!(bnd::Boundary, dis::Discretization{Polar}, sys::System, ε, F, r; display::Bool=false)
    R = broadcast((a,b)->a,dis.x...)
    Θ = broadcast((a,b)->b,dis.x...)
    sub_pixel_num = dis.sub_pixel_num
    if dis.N[2] == 1

    else
        sub_x = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        sub_y = Array{Float64}(undef, sub_pixel_num, sub_pixel_num)
        domains = Array{Int}(undef, sub_pixel_num, sub_pixel_num)
        ε_temp = Array{ComplexF64}(undef,sub_pixel_num, sub_pixel_num)
        F_temp = Array{Float64}(undef,sub_pixel_num, sub_pixel_num)
        if display
            pg = Progress((size(x,1)-2)*(size(y,2)-2), PROGRESS_UPDATE_TIME::Float64, "sub-pixel smoothing ")
        end
        for h ∈ CartesianIndices(r)
            i,j = h[1],h[2]
            if all((1,1).<(i,j).<size(r))
                nearestNeighborFlag = r[i,j]!==r[i,j+1] || r[i,j]!==r[i,j-1] || r[i,j]!==r[i+1,j] || r[i,j]!==r[i-1,j]
                nextNearestNeighborFlag = r[i,j]!==r[i+1,j+1] || r[i,j]!==r[i-1,j-1] || r[i,j]!==r[i+1,j-1] || r[i,j]!==r[i-1,j+1]
                if nearestNeighborFlag || nextNearestNeighborFlag
                    r_min = (R[i,j]+R[i-1,j])/2
                    r_max = (R[i,j]+R[i+1,j])/2
                    θ_min = (Θ[i,j]+Θ[i,j-1])/2
                    θ_max = (Θ[i,j]+Θ[i,j+1])/2
                    for kj ∈ 1:sub_pixel_num
                        sub_θ = θ_min + (kj-1)*(θ_max-θ_min)/(sub_pixel_num-1)
                        sc = sincos(sub_θ)
                        for ki ∈ 1:sub_pixel_num
                            sub_r = r_min + (ki-1)*(r_max-r_min)/(sub_pixel_num-1)
                            sub_x[ki,kj] = sub_r*sc[2]
                            sub_y[ki,kj] = sub_r*sc[1]
                        end
                    end
                    which_domains!(domains, sub_x, sub_y, bnd, sys)
                    for l ∈ eachindex(domains)
                        d = domains[l]
                        ε_temp[l] = sys.ε_by_region[d](sub_x[d],sub_y[d],sys.params_by_region[d])
                        F_temp[l] = sys.F_by_region[d](sub_x[d],sub_y[d],sys.params_by_region[d])
                    end
                    ε[i,j] = mean(ε_temp)
                    F[i,j] = mean(F_temp)
                end
                if display
                    next!(pg)
                end
            end
        end
    end
    return nothing
end


function which_domains(x,y,bnd,sys)
    domains = zeros(Int,size(x)...)
    which_domains!(domains,x,y,bnd,sys)
    return domains
end
function which_domains!(domains,x,y,bnd,sys) # need to add boundary
    for i ∈ eachindex(x)
        domains[i] = which_domain(x[i],y[i],bnd,sys)
    end
    return nothing
end


"""
    which_domain(x, y, bnd, sys)
"""
function which_domain(x, y, bnd::Boundary, sys::System)
    domain = findfirst(map(z->z.is_in_domain(x,y),sys.domains))
    @assert !isnothing(domain) "domain not found for (x,y)=($x,$y). has a background domain been defined?"
    return domain
end


################################################################################
###############  STANDARD DOMAINS  #######################################
################################################################################
"""
    left_domain(x, y, domain_index, bnd, sys)
"""
left_domain(x,y,idx::Int,bnd::Boundary,sys::System) =
    sys.domains[idx].which_asymptote == :left && x<bnd.∂Ω_tr[1][1]


"""
    right_domain(x, y, domain_index, bnd, sys)
"""
function right_domain(x,y,idx::Int,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :right && x>bnd.∂Ω_tr[1][2]
end


"""
    bottom_domain(x, y, domain_index, bnd, sys)
"""
function bottom_domain(x,y,idx::Int,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :bottom && y<bnd.∂Ω_tr[2][1]
end


"""
    top_domain(x, y, domain_index, bnd, sys)
"""
function top_domain(x,y,idx::Int,bnd::Boundary,sys::System)
    return sys.domains[idx].which_asymptote == :top && y>bnd.∂Ω_tr[2][2]
end
