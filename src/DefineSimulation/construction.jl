# TODO: sub_pixel_smoothing in 2-dim case, take care of smoothing along edges
# TODO: fix Bravais in 1-dimension

################################################################################
###############   SUB PIXEL SMOOTHING    #######################################
################################################################################
function construct_εFr(bnd::Boundary, dis::Discretization, sys::System)
    x = dis.X
    y = dis.Y
    domains = Array{Int}(undef,size(x)...)
    r = domains
    which_domains!(domains, x, y, bnd, sys)

    xb, yb = Array{Float64}(undef,size(x)...), Array{Float64}(undef,size(y)...)
    bravais_coordinates_unit_cell!(xb,yb,x,y,domains,sys)

    n1, n2, f = [1.0], [0.0], [0.0]
    ε = Array{ComplexF64}(undef,dis.N[1],dis.N[2])
    F = Array{Float64}(undef,dis.N[1],dis.N[2])
    for i ∈ eachindex(r)
        j = r[i]
        ε[i] = sys.ε_by_region[j](x[j],y[j],sys.params_by_region[r[j]],n1,n2)
        F[i] = sys.F_by_region[j](x[j],y[j],sys.params_by_region[r[j]],f)
    end
    return ε, F, r
end


"""
    ε, regions = sub_pixel_smoothing(bnd, dis, sys; display=false)
"""
function sub_pixel_smoothing!(bnd::Boundary, dis::Discretization, sys::System, ε, F, r; display::Bool=false)
    x = dis.X
    y = dis.Y
    sub_pixel_num = dis.sub_pixel_num
    n1, n2, f = [1.0], [0.0], [0.0]
    if dis.N[2] == 1
        # sub_x = Array{Float64}(undef, sub_pixel_num,1)
        # sub_y = [y[1]]
        # xb = Array{Float64}(undef, sub_pixel_num, 1)
        # yb = Array{Float64}(undef, sub_pixel_num, 1)
        # sub_regions = Array{Int}(undef, sub_pixel_num, 1)
        # sub_domains = Array{Int}(undef, sub_pixel_num, 1)
        # for i ∈ 2:(size(x,1)-1)
        #     nearestNeighborFlag = r[i]!==r[i+1] || r[i]!==r[i-1]
        #     if nearestNeighborFlag
        #         x_min = (x[i]+x[i-1])/2
        #         x_max = (x[i]+x[i+1])/2
        #         sub_x[:] = LinRange(x_min, x_max, sub_pixel_num)
        #         sub_domains[:] = which_domain.(sub_x, sub_y, Ref(bnd), Ref(sys))
        #         bravais_coordinates_unit_cell!(xb, yb, sub_x, sub_y, sub_domains, sys)
        #         sub_regions[:] = which_region.(xb, yb, sub_domains, Ref(sys))
        #         ε[i] = mean([sys.ε_by_region[sub_regions[i]](sub_x[i],sub_y[1],sys.params_by_region[sub_regions[i]],n1,n2) for i ∈ eachindex(sub_regions)])
        #         F[i] = mean([sys.F_by_region[sub_regions[i]](sub_x[i],sub_y[1],sys.params_by_region[sub_regions[i]],f) for i ∈ eachindex(sub_regions)])
        #     end
        # end
    elseif dis.N[1] == 1
        # sub_x = [x[1]]
        # sub_y = Array{Float64}(undef, 1, sub_pixel_num)
        # xb = Array{Float64}(undef, 1, sub_pixel_num)
        # yb = Array{Float64}(undef, 1, sub_pixel_num)
        # sub_regions = Array{Int}(undef, 1, sub_pixel_num)
        # sub_domains = Array{Int}(undef, 1, sub_pixel_num)
        # for i ∈ 2:(size(y,2)-1)
        #     nearestNeighborFlag = r[i]!==r[i+1] || r[i]!==r[i-1]
        #     if nearestNeighborFlag
        #         y_min = (y[i]+y[i-1])/2
        #         y_max = (y[i]+y[i+1])/2
        #         sub_y[:] = LinRange(y_min, y_max, sub_pixel_num)
        #         sub_domains[:] = which_domain.(sub_x, sub_y, Ref(bnd), Ref(sys))
        #         bravais_coordinates_unit_cell!(xb, yb, sub_x, sub_y, sub_domains, sys)
        #         sub_regions[:] = which_region.(xb, yb, sub_domains, Ref(sys))
        #         ε[i] = mean([sys.ε_by_region[sub_regions[i]](sub_x[1],sub_y[i],sys.params_by_region[sub_regions[i]],n1,n2) for i ∈ eachindex(sub_regions)])
        #         F[i] = mean([sys.F_by_region[sub_regions[i]](sub_x[1],sub_y[i],sys.params_by_region[sub_regions[i]],f) for i ∈ eachindex(sub_regions)])
        #     end
        # end
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
        for i ∈ 2:(size(x,1)-1), j ∈ 2:(size(y,2)-1)
            nearestNeighborFlag = r[i,j]!==r[i,j+1] || r[i,j]!==r[i,j-1] || r[i,j]!==r[i+1,j] || r[i,j]!==r[i-1,j]
            nextNearestNeighborFlag = r[i,j]!==r[i+1,j+1] || r[i,j]!==r[i-1,j-1] || r[i,j]!==r[i+1,j-1] || r[i,j]!==r[i-1,j+1]
            if nearestNeighborFlag || nextNearestNeighborFlag
                x_min = (x[i]+x[i-1])/2
                y_min = (y[j]+y[j-1])/2
                x_max = (x[i]+x[i+1])/2
                y_max = (y[j]+y[j+1])/2
                for k ∈ CartesianIndices(sub_x)
                    sub_x[k] = x_min + (k[1]-1)*(x_max-x_min)/(sub_pixel_num-1)
                    sub_y[k] = y_min + (k[2]-1)*(y_max-y_min)/(sub_pixel_num-1)
                end
                which_domains!(domains, sub_x, sub_y, bnd, sys)
                bravais_coordinates_unit_cell!(xb, yb, sub_x, sub_y, domains, sys)
                for i ∈ eachindex(xb)
                    j = domains[i]
                    ε_temp[i] = sys.ε_by_region[j](sub_x[j],sub_y[j],sys.params_by_region[j],n1,n2)
                    F_temp[i] = sys.F_by_region[j](sub_x[j],sub_y[j],sys.params_by_region[j],f)
                end
                ε[i,j] = mean(ε_temp)
                F[i,j] = mean(F_temp)
            end
            if display
                next!(pg)
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
