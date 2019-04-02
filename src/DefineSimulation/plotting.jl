struct PaintBucket{Tbnd}
    N::Int64
    x_f_start::Float64
    x_f_stop::Float64
    x_f::LinRange{Float64}
    x_t_start::Float64
    x_t_stop::Float64
    x_t::LinRange{Float64}
    y_f_start::Float64
    y_f_stop::Float64
    y_f::LinRange{Float64}
    y_t_start::Float64
    y_t_stop::Float64
    y_t::LinRange{Float64}
    ε::Array{Complex{Float64},2}
    F::Array{Float64,2}
    ψ_plot::Array{Complex{Float64},3}
    boundary_data::Tbnd
    pml_data::NTuple{4,Tuple{Array{Float64,1},Array{Float64,1}}}

    function PaintBucket(sim::Simulation{Tsys,Tbnd,Discretization{Polar}}, ψ, truncate::Bool) where {Tsys,Tbnd}

        r_f = sim.dis.x[1][1] .+ sim.dis.dx[1]*(-1:sim.dis.N[1])
        θ_f = sim.dis.x[2][1] .+ sim.dis.dx[2]*(-1:sim.dis.N[2])

        x_f_start = minimum(sim.dis.X)-sim.dis.dx[1]/2
        x_f_stop  = maximum(sim.dis.X)+sim.dis.dx[1]/2
        x_f = LinRange(x_f_start, x_f_stop, 2sim.dis.N[1])

        y_f_start = minimum(sim.dis.Y)-sim.dis.dx[1]/2
        y_f_stop  = maximum(sim.dis.Y)+sim.dis.dx[1]/2
        y_f = LinRange(y_f_start, y_f_stop, 2sim.dis.N[1])

        ε_temp = Array{ComplexF64}(undef, sim.dis.N[1]+2, sim.dis.N[2]+2)
        ε_temp[2:end-1,2:end-1] = sim.sys.ε
        @views ε_temp[:,1] = ε_temp[:,end-1]
        @views ε_temp[:,end] = ε_temp[:,2]
        @views ε_temp[1,:] = ε_temp[2,:]
        @views ε_temp[end,:] = ε_temp[end-1,:]
        # ε_itp = extrapolate(scale(interpolate(ε_temp,BSpline(Quadratic(Periodic(OnGrid())))),r_f,θ_f),complex(NaN,NaN))
        ε_itp = extrapolate(scale(interpolate(ε_temp,BSpline(Linear())),r_f,θ_f),complex(NaN,NaN))
        ε = @. ε_itp(hypot(x_f,y_f'),sim.dis.x[2][1]+mod(atan(y_f',x_f)-sim.dis.x[2][1],2π))

        F_temp = Array{Float64}(undef, sim.dis.N[1]+2, sim.dis.N[2]+2)
        F_temp[2:end-1,2:end-1] = sim.sys.F
        @views F_temp[:,1] = F_temp[:,end-1]
        @views F_temp[:,end] = F_temp[:,2]
        @views F_temp[1,:] = F_temp[2,:]
        @views F_temp[end,:] = F_temp[end-1,:]
        # F_itp = extrapolate(scale(interpolate(F_temp,BSpline(Cubic(Periodic(OnGrid())))),r_f,θ_f),NaN)
        F_itp = extrapolate(scale(interpolate(F_temp,BSpline(Linear())),r_f,θ_f),NaN)
        F = @. F_itp(hypot(x_f,y_f'),sim.dis.x[2][1]+mod(atan(y_f',x_f)-sim.dis.x[2][1],2π))

        N, M, X, Y, x_t, y_t, ∂Ω, ψ = apply_truncation(sim, ψ, truncate)

        r_t = x_t[1] .+ sim.dis.dx[1]*(-1:M[1])
        θ_t = y_t[1] .+ sim.dis.dx[2]*(-1:M[2])

        x_t_start = minimum(X)-sim.dis.dx[1]/2
        x_t_stop  = maximum(X)+sim.dis.dx[1]/2
        x_t = LinRange(x_t_start, x_t_stop, 2M[1])

        y_t_start = minimum(Y)-sim.dis.dx[1]/2
        y_t_stop  = maximum(Y)+sim.dis.dx[1]/2
        y_t = LinRange(y_t_start, y_t_stop, 2M[1])

        # root_r = sqrt.(hypot.(X,Y))
        ψ_temp = Array{ComplexF64}(undef, M[1]+2, M[2]+2)
        ψ_plot = Array{ComplexF64}(undef, 2M[1], 2M[1], N)
        for i ∈ 1:N
            @views ψ_temp[2:end-1,2:end-1] = reshape(ψ[:,i], M[1], M[2])#./root_r
            @views ψ_temp[:,1] = ψ_temp[:,end-1]
            @views ψ_temp[:,end] = ψ_temp[:,2]
            @views ψ_temp[1,:] = ψ_temp[2,:]
            @views ψ_temp[end,:] = ψ_temp[end-1,:]
            ψ_itp = extrapolate(scale(interpolate(ψ_temp,BSpline(Cubic(Periodic(OnGrid())))),r_t,θ_t),complex(NaN,NaN))
            ψ_plot[:,:,i] = @. ψ_itp(hypot(x_t,y_t'),sim.dis.x[2][1]+mod(atan(y_t',x_t)-sim.dis.x[2][1],2π))
        end

        boundary_data1_f = (r_f[1]*cos.(θ_f),r_f[1]*sin.(θ_f))
        boundary_data2_f = (r_f[end]*cos.(θ_f),r_f[end]*sin.(θ_f))
        boundary_data3_f = (r_f*cos(θ_f[1]),r_f*sin(θ_f[1]))
        boundary_data4_f = (r_f*cos(θ_f[end]),r_f*sin(θ_f[end]))

        N_pml = 101
        pml_r_1 = vcat(
                    LinRange(sim.dis.x[1][1], sim.dis.x_tr[1][1], N_pml),
                    fill(sim.dis.x_tr[1][1], N_pml),
                    LinRange(sim.dis.x_tr[1][1], sim.dis.x[1][1], N_pml),
                fill(sim.dis.x[1][1],N_pml)
                )
        pml_θ_1 = vcat(
                    fill(sim.dis.x[2][1],N_pml),
                    LinRange(sim.dis.x[2][1], sim.dis.x[2][end], N_pml),
                    fill(sim.dis.x[2][end],N_pml),
                    LinRange(sim.dis.x[2][end], sim.dis.x[2][1], N_pml)
                    )

        pml_r_2 = vcat(
                    LinRange(sim.dis.x_tr[1][end], sim.dis.x[1][end], N_pml),
                    fill(sim.dis.x[1][end], N_pml),
                    LinRange(sim.dis.x[1][end], sim.dis.x_tr[1][end], N_pml),
                    fill(sim.dis.x_tr[1][end],N_pml)
                    )
        pml_θ_2 = vcat(
                    fill(sim.dis.x[2][1],N_pml),
                    LinRange(sim.dis.x[2][1], sim.dis.x[2][end], N_pml),
                    fill(sim.dis.x[2][end],N_pml),
                    LinRange(sim.dis.x[2][end], sim.dis.x[2][1], N_pml)
                    )

        pml_r_3 = vcat(
                    LinRange(sim.dis.x[1][1], sim.dis.x[1][end], N_pml),
                    fill(sim.dis.x[1][end], N_pml),
                    LinRange(sim.dis.x[1][end], sim.dis.x[1][1], N_pml),
                    fill(sim.dis.x[1][1],N_pml)
                    )
        pml_θ_3 = vcat(
                    fill(sim.dis.x[2][1],N_pml),
                    LinRange(sim.dis.x[2][1], sim.dis.x_tr[2][1], N_pml),
                    fill(sim.dis.x_tr[2][1],N_pml),
                    LinRange(sim.dis.x_tr[2][1], sim.dis.x[2][1], N_pml)
                    )

        pml_r_4 = vcat(
                    LinRange(sim.dis.x[1][1], sim.dis.x[1][end], N_pml),
                    fill(sim.dis.x[1][end], N_pml),
                    LinRange(sim.dis.x[1][end], sim.dis.x[1][1], N_pml),
                    fill(sim.dis.x[1][1],N_pml)
                    )
        pml_θ_4 = vcat(
                    fill(sim.dis.x_tr[2][end],N_pml),
                    LinRange(sim.dis.x_tr[2][end], sim.dis.x[2][end], N_pml),
                    fill(sim.dis.x[2][end],N_pml),
                    LinRange(sim.dis.x[2][end], sim.dis.x_tr[2][end], N_pml)
                    )

        pml_data1 = (pml_r_1.*cos.(pml_θ_1),pml_r_1.*sin.(pml_θ_1))
        pml_data2 = (pml_r_2.*cos.(pml_θ_2),pml_r_2.*sin.(pml_θ_2))
        pml_data3 = (pml_r_3.*cos.(pml_θ_3),pml_r_3.*sin.(pml_θ_3))
        pml_data4 = (pml_r_4.*cos.(pml_θ_4),pml_r_4.*sin.(pml_θ_4))

        bnd_data = (boundary_data1_f, boundary_data2_f, boundary_data3_f, boundary_data4_f)

        return new{typeof(bnd_data)}(N, x_f_start, x_f_stop, x_f, x_t_start, x_t_stop, x_t,
            y_f_start, y_f_stop, y_f, y_t_start, y_t_stop, y_t,
            ε, F, ψ_plot,
            bnd_data,
            (pml_data1, pml_data2, pml_data3, pml_data4) )
    end


    function PaintBucket(sim::Simulation{Tsys,Tbnd,Discretization{Cartesian}}, ψ, truncate::Bool) where Tsys<:System where Tbnd<:Boundary

        N, M, X, Y, x_t, y_t, ∂Ω, ψ = apply_truncation(sim, ψ, truncate)

        x_f_start = sim.bnd.∂Ω[1][1]
        x_f_stop  = sim.bnd.∂Ω[1][2]
        x_f = sim.dis.x[1][:]

        y_f_start = sim.bnd.∂Ω[2][1]
        y_f_stop  = sim.bnd.∂Ω[2][2]
        y_f = sim.dis.x[2][:]

        x_t_start = ∂Ω[1][1]
        x_t_stop  = ∂Ω[1][2]

        y_t_start = ∂Ω[2][1]
        y_t_stop  = ∂Ω[2][2]

        ε = sim.sys.ε
        # ε = (1 .+ 1im*sim.sys.Σe/10).*sim.sys.ε
        F = sim.sys.F

        ψ_plot = Array{ComplexF64}(undef, M[1], M[2], N)
        for i in 1:N
            ψ_plot[:,:,i] = reshape(ψ[:,i], M[1], M[2])
        end

        boundary_data1_f = ([x_f[1],x_f[1]],[y_f[1],y_f[end]])
        boundary_data2_f = ([x_f[end],x_f[end]],[y_f[1],y_f[end]])
        boundary_data3_f = ([x_f[1],x_f[end]],[y_f[1],y_f[1]])
        boundary_data4_f = ([x_f[1],x_f[end]],[y_f[end],y_f[end]])

        N_pml = 2
        pml_x_1 = vcat(
                    LinRange(sim.dis.x[1][1], sim.dis.x_tr[1][1], N_pml),
                    fill(sim.dis.x_tr[1][1], N_pml),
                    LinRange(sim.dis.x_tr[1][1], sim.dis.x[1][1], N_pml),
                    fill(sim.dis.x[1][1],N_pml)
                    )
        pml_y_1 = vcat(
                    fill(sim.dis.x[2][1],N_pml),
                    LinRange(sim.dis.x[2][1], sim.dis.x[2][end], N_pml),
                    fill(sim.dis.x[2][end],N_pml),
                    LinRange(sim.dis.x[2][end], sim.dis.x[2][1], N_pml)
                    )

        pml_x_2 = vcat(
                    LinRange(sim.dis.x_tr[1][end], sim.dis.x[1][end], N_pml),
                    fill(sim.dis.x[1][end], N_pml),
                    LinRange(sim.dis.x[1][end], sim.dis.x_tr[1][end], N_pml),
                    fill(sim.dis.x_tr[1][end],N_pml)
                    )
        pml_y_2 = vcat(
                    fill(sim.dis.x[2][1],N_pml),
                    LinRange(sim.dis.x[2][1], sim.dis.x[2][end], N_pml),
                    fill(sim.dis.x[2][end],N_pml),
                    LinRange(sim.dis.x[2][end], sim.dis.x[2][1], N_pml)
                    )

        pml_x_3 = vcat(
                    LinRange(sim.dis.x[1][1], sim.dis.x[1][end], N_pml),
                    fill(sim.dis.x[1][end], N_pml),
                    LinRange(sim.dis.x[1][end], sim.dis.x[1][1], N_pml),
                    fill(sim.dis.x[1][1],N_pml)
                    )
        pml_y_3 = vcat(
                    fill(sim.dis.x[2][1],N_pml),
                    LinRange(sim.dis.x[2][1], sim.dis.x_tr[2][1], N_pml),
                    fill(sim.dis.x_tr[2][1],N_pml),
                    LinRange(sim.dis.x_tr[2][1], sim.dis.x[2][1], N_pml)
                    )

        pml_x_4 = vcat(
                    LinRange(sim.dis.x[1][1], sim.dis.x[1][end], N_pml),
                    fill(sim.dis.x[1][end], N_pml),
                    LinRange(sim.dis.x[1][end], sim.dis.x[1][1], N_pml),
                    fill(sim.dis.x[1][1],N_pml)
                    )
        pml_y_4 = vcat(
                    fill(sim.dis.x_tr[2][end],N_pml),
                    LinRange(sim.dis.x_tr[2][end], sim.dis.x[2][end], N_pml),
                    fill(sim.dis.x[2][end],N_pml),
                    LinRange(sim.dis.x[2][end], sim.dis.x_tr[2][end], N_pml)
                    )

        pml_data1 = (pml_x_1,pml_y_1)
        pml_data2 = (pml_x_2,pml_y_2)
        pml_data3 = (pml_x_3,pml_y_3)
        pml_data4 = (pml_x_4,pml_y_4)

        M[1]==1 ? x_t = y_t : nothing
        M[1]==1 ? x_f_start = y_f_start : nothing
        M[1]==1 ? x_f_stop = y_f_stop : nothing
        M[1]==1 ? x_f = y_f : nothing
        M[1]==1 ? x_t_start = y_t_start : nothing
        M[1]==1 ? x_t_stop = y_t_stop : nothing

        bnd_data = (boundary_data1_f, boundary_data2_f, boundary_data3_f, boundary_data4_f)
        return new{typeof(bnd_data)}(N, x_f_start, x_f_stop, x_f, x_t_start, x_t_stop, x_t[:],
            y_f_start, y_f_stop, y_f, y_t_start, y_t_stop, y_t[:],
            ε, F, ψ_plot,
            bnd_data,
            (pml_data1, pml_data2, pml_data3, pml_data4) )
    end
end


"""
    p = plot(domain)
"""
@recipe function f(dom::Domain)
    (dom.is_in_domain)
end


"""
    p = plot(sim; by=nothing, truncate=false)
"""
@recipe function f(sim::Simulation; by=nothing, truncate=false)
    pb = PaintBucket(sim, Array{ComplexF64}(undef,prod(sim.dis.N),1), truncate)
    if isnothing(by)
        layout --> (1,3)
        (sim,pb,by,1)
    else
        (sim,pb,by,1,1)
    end
end


# internal use, sim summary
@recipe function f(sim::Simulation, pb::PaintBucket, by, _1)
    bys = [:real, :imag, :F]
    for i ∈ 1:3
        @series begin
            subplot := i
            (sim, pb, bys[i], 1, 1)
        end
    end
end


# interal use only, sim elementary
@recipe function f(sim::Simulation, pb::PaintBucket, by, _1, _2)
    aspect_ratio --> 1
    grid --> false
    xlims --> [pb.x_t_start,pb.x_t_stop]
    ylims --> [pb.y_t_start,pb.y_t_stop]
    colorbar --> false
    overwrite_figure --> false
    levels --> 15
    lw --> 6
    seriestype --> :heatmap
    legend --> false
    framestyle --> :grid

    if by ∈ [:real, real, "real"]
        color --> :sequential
        title --> "Re n(x)"
        p = real(sqrt.(pb.ɛ))
        pm = p[.!isnan.(p)]
        clims := (minimum(pm), maximum(pm))
    elseif by ∈ [:imag, imag, "imag"]
        color --> :diverging
        title --> "Im n(x)"
        p = imag(sqrt.(pb.ɛ))
        pm = p[.!isnan.(p)]
        clims := (-maximum(abs.(pm)), maximum(abs.(pm)))
    elseif by ∈ [:F,"F"]
        color --> :diverging
        title --> "F(x)"
        p = pb.F
        pm = p[.!isnan.(p)]
        clims := (-maximum(abs.(pm)), maximum(abs.(pm)))
    else
        throw(ArgumentError("unrecognized `by` keyword"))
    end
    @series begin
        pb.x_f, pb.y_f, permutedims(p)
    end
    for k ∈ 1:2 # cycle over dims
        for i ∈ eachindex(sim.bnd.bc) # cycle over sides
            if isDirichlet(sim.bnd.bc[k][i])
                linestyle := :solid
            elseif isNeumann(sim.bnd.bc[k][i])
                linestyle := :dash
            end
            if !isnoBL(sim.bnd.bl[k][i])
                if isPML(sim.bnd.bl[k][i])
                    fill_color = :red
                elseif iscPML(sim.bnd.bl[k][i])
                    fill_color = :blue
                end
                @series begin
                    lw := 0
                    seriestype := :path
                    color := fill_color
                    alpha := .25
                    fill --> (0,.5,fill_color)
                    pb.pml_data[i]
                end
            end
            if isDirichlet(sim.bnd.bc[k][i]) | isNeumann(sim.bnd.bc[k][i])
                @series begin
                    lw := 2
                    color := :black
                    seriestype := :path
                    pb.boundary_data[i]
                end
            end
        end
    end
end


"""
    p = plot(sim, ψ, [inds; by=nothing, summary=false, truncate=false])
"""
@recipe function f(sim::Simulation, ψ::AbstractArray, inds=1:size(ψ,2); by=nothing, structure=false, truncate=false)
    pb = PaintBucket(sim, ψ[:,inds], truncate)
    n = length(inds)
    if isnothing(by)
        if structure
            layout --> (1+n,3)
        else
            layout --> (n,3)
        end
    else
        if structure
            layout --> (3+n)
        else
            layout --> (n)
        end
    end
    if structure
        @series (sim,pb,by,1)
    end
    if isnothing(by)
        for i ∈ 1:n
            @series begin
                subplot := 3*(structure+i-1)+1
                (sim,pb,i,:real,1,1,1)
            end
            @series begin
                subplot := 3*(structure+i-1)+2
                (sim,pb,i,:imag,1,1,1)
            end
            @series begin
                subplot := 3*(structure+i-1)+3
                (sim,pb,i,:abs2,1,1,1)
            end
        end
    else
        for i ∈ 1:n
            @series begin
                subplot := 3*structure+i
                (sim,pb,i,by,1,1,1)
            end
        end
    end
end


# elementary plot sim of ψ[:,:,i]
@recipe function f(sim::Simulation,pb::PaintBucket,i::Int,by,_1,_2,_3)
    aspect_ratio --> 1
    grid --> false
    xlims --> [pb.x_t_start,pb.x_t_stop]
    ylims --> [pb.y_t_start,pb.y_t_stop]
    colorbar --> false
    overwrite_figure --> false
    levels --> 15
    lw --> 6
    seriestype --> :heatmap
    legend --> false
    framestyle --> :grid

    p = pb.ψ_plot[:,:,i]
    if by ∈ [:real,real,"real"]
        color --> :diverging
        title --> "Re psi(x)"
        p = real(p)
    elseif by ∈ [:imag,imag,"imag"]
        color --> :diverging
        title --> "Im psi(x)"
        p = imag(p)
    elseif by ∈ [:abs,abs,"abs"]
        color --> :sequential
        title --> "|psi(x)|"
        p = abs.(p)
    elseif by ∈ [:abs2,abs2,"abs2"]
        color --> :sequential
        title --> "|psi(x)|²"
        p = abs2.(p)
    else
        throw(ArgumentError("unrecognized `by` keyword"))
    end
    @series begin
        pb.x_t, pb.y_t, permutedims(p)
    end
end


function apply_truncation(sim::Simulation, ψ, truncate::Bool)
    if isempty(ψ)
        N=0
    else
        N = size(ψ,2)
    end
    if !isempty(ψ)
        idx = truncate ? sim.dis.X_idx : 1:prod(sim.dis.N)
    else
        idx = Int[]
    end
    if !truncate || iszero(N)
        M = sim.dis.N
        X, Y = sim.dis.X, sim.dis.Y
        x, y = sim.dis.x[1], sim.dis.x[2]
        ∂Ω = sim.bnd.∂Ω
    else
        M = sim.dis.N_tr
        X = reshape(sim.dis.X[sim.dis.X_idx],M[1],M[2])
        Y = reshape(sim.dis.Y[sim.dis.X_idx],M[1],M[2])
        x, y = sim.dis.x_tr[1], sim.dis.x_tr[2]
        ∂Ω = sim.bnd.∂Ω_tr
        ψ = ψ[idx,:]
    end
    return N, M, X, Y, x, y, ∂Ω, ψ
end


# ################################################################################
# ########## ANIMATION
# ################################################################################
# """
#     iterator = wave(sim, ψ; by=real, n=60, seriestype=:heatmap)
#
# input for Plots.animate
#
# Use cases:
#
# `animate(wave(sim,ψ), file_name)` creates a .gif with filename
#
# `animate(wave(sim,ψ; n=20), file_name, fps=10)` creates a 2 second movie
#
# Note: default `fps`=20, and `n`=60, so default movie is 3 seconds long
# """
# function wave(sim::Simulation, ψ; truncate=true, by=real, n=60)
#     if 1 ∈ sim.dis.N
#         return imap( ϕ->(sim, exp(-1im*ϕ)*ψ, by, truncate, 1), 0:2π/n:2π*(1-1/n))
#     else
#         return imap( ϕ->(sim, exp(-1im*ϕ)*ψ, by, truncate, 1, 1), 0:2π/n:2π*(1-1/n))
#     end
# end
