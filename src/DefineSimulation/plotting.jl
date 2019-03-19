#TODO: check 1d plots, to waveguixe_dispersion

################################################################################
########## SIMULATION
################################################################################
"""
    p = plot(sim)
"""
@recipe function f(sim::Simulation)
    (sim, ComplexF64[])
end


################################################################################
########## SOLUTIONS
################################################################################
"""
    p = plot(sim, ψ; by=nothing, truncate=true)

plots
to turn off translucent effect, add optional argument `seriesalpha=0`
vary type of plot with `seriestype`, e.g. `seriestype=:contour`
"""
@recipe function f(sim::Simulation, ψ::AbstractArray; by=nothing, truncate=true)
    if 1 ∈ sim.dis.N
        (sim, ψ, by, truncate, 1)
    else
        (sim, ψ, by, truncate, 2, 2)
    end
end


# 1d plot
@recipe function f(sim::Simulation, ψ::AbstractArray, by::Union{Function,Nothing}, truncate::Bool, dim1::Int)

    (N, x_f_start, x_f_stop, x_f, x_t_start, x_t_stop, x_t,
        y_f_start, y_f_stop, y_f, y_t_start, y_t_stop, y_t,
        ε, F, ψ_plot, boundary_data, pml_data) = prep_data(sim, ψ, truncate)

    cmapc, cmapk, cmapsim1, cmapsim2, n_mult, F_sign = fix_colormap(COLOR_SCHEME)

    xlims --> [x_t_start,x_t_stop]
    colorbar --> false
    overwrite_figure --> false
    levels --> 15
    lw --> 2
    seriestype := :line
    legend --> false

    if sim.dis.N[1]==1
        bottom_label = "y"
    else
        bottom_label = "x"
    end

    if by==nothing

        layout --> (1+N,3)

        # plot simulation
        titles = ["real", "imag", "F/abs²"]
        @series begin
            subplot := 1
            title := titles[1]
            seriestype := :line
            n₁ = real(sqrt.(ɛ-1im*sim.tls.D₀*F))
            x_f[:], n₁[:]
        end
        @series begin
            subplot := 2
            title := titles[2]
            seriestype := :line
            n₂ = imag(sqrt.(ɛ-1im*sim.tls.D₀*F))
            x_f[:], n₂[:]
        end
        @series begin
            subplot := 3
            title := titles[3]
            seriestype := :line
            x_f[:], F[:]
        end

        for i ∈ 1:N
            ψ = ψ_plot[:,:,i]
            @series begin
                subplot := 3i+1
                x_t[:], real(ψ[:])
            end
            @series begin
                subplot := 3i+2
                x_t[:], imag(ψ[:])
            end
            @series begin
                subplot := 3i+3
                x_t[:], abs2.(ψ[:])
            end
        end
    else
        if iszero(N)
            layout --> (1,1)
            @series begin
                subplot := 1
                n = sqrt.(ɛ-1im*sim.tls.D₀*F)
                x_f[:], by.(n[:])
            end
        else
            if round(Int,N/3)==N/3
                N_col=3
                N_row = ceil(Int,N/3)
            elseif round(Int,N/2)==N/2
                N_col=2
                N_row = ceil(Int,N/2)
            elseif N==1
                N_col=1
                N_row=1
            else
                N_col = 3
                N_row = ceil(Int,N/3)
            end
            layout --> (N_row,N_col)
            for i ∈ 1:N
                ψ = ψ_plot[:,:,i]
                @series begin
                    subplot := i
                    x_t[:], by.(ψ[:])
                end
            end
        end
    end

    # if by==nothing
    #
    #     size --> 250*[3,1+N]
    #     layout := (1+N,3)
    #
    #     @series begin
    #         ylabel := string("Re{n(", bottom_label, ")}")
    #         xlabel --> bottom_label
    #         subplot := 1
    #         n = real(sqrt.(ɛ-1im*sim.tls.D₀*F))
    #         ylims := (minimum(n), n_mult*maximum(n))
    #         x, n
    #     end
    #     @series begin
    #         ylabel := string("Im{n(", bottom_label, ")}")
    #         xlabel --> bottom_label
    #         subplot := 2
    #         n = imag(sqrt.(ɛ-1im*sim.tls.D₀*F))
    #         ylims := (-maximum(abs.(n)), maximum(abs.(n)))
    #         x, n
    #     end
    #     @series begin
    #         ylabel := string("F(", bottom_label, ")")
    #         xlabel --> bottom_label
    #         subplot := 3
    #         ylims := (-maximum(abs.(F)), maximum(abs.(F)))
    #         x, F
    #     end
    #
    #     for i ∈ 1:N
    #         ψ = ψ_plot[:,i]
    #         @series begin
    #             ylabel := LaTeXString("Re\\{\\psi\\}")
    #             xlabel --> bottom_label
    #             subplot := 3i+1
    #             ylims = (-maximum(abs.(ψ)), maximum(abs.(ψ)))
    #             x, real(ψ)
    #         end
    #         @series begin
    #             ylabel := LaTeXString("Im\\{\\psi\\}")
    #             xlabel --> bottom_label
    #             subplot := 3i+2
    #             ylims = (-maximum(abs.(ψ)), maximum(abs.(ψ)))
    #             x, imag(ψ)
    #         end
    #         @series begin
    #             ylabel := LaTeXString("|\\psi|^2")
    #             xlabel --> bottom_label
    #             subplot := 3i+3
    #             ylims := (0, maximum(abs2.(ψ)))
    #             x, abs2.(ψ)
    #         end
    #     end
    # else
    #     layout := N
    #     for i ∈ 1:N
    #         ψ = ψ_plot[:,i]
    #         @series begin
    #             if by ∈ [abs, abs2]
    #                 ylims := (0, by.(ψ))
    #             else
    #                 ylims := (-maximum(abs.(ψ)), maximum(abs.(ψ)))
    #             end
    #             x, by.(ψ)
    #         end
    #     end
    # end
end


# 2d plot
@recipe function f(sim::Simulation, ψ::AbstractArray, by::Union{Function,Nothing}, truncate::Bool, dim1::Int, dim2::Int)

    (N, x_f_start, x_f_stop, x_f, x_t_start, x_t_stop, x_t,
        y_f_start, y_f_stop, y_f, y_t_start, y_t_stop, y_t,
        ε, F, ψ_plot, boundary_data, pml_data) = prep_data(sim, ψ, truncate)

    cmapc, cmapk, cmapsim1, cmapsim2, n_mult, F_sign = fix_colormap(COLOR_SCHEME)

    aspect_ratio --> 1
    xlims --> [x_t_start,x_t_stop]
    ylims --> [y_t_start,y_t_stop]
    colorbar --> false
    overwrite_figure --> false
    levels --> 15
    lw --> 6
    seriestype --> :heatmap
    legend --> false

    if by==nothing

        layout --> (1+N,3)

        # plot simulation
        titles = ["real", "imag", "F/abs²"]
        @series begin
            subplot := 1
            title := titles[1]
            color := cmapsim1
            seriestype := :heatmap
            n₁ = real(sqrt.(ɛ-1im*sim.tls.D₀*F))
            clims := (minimum(n₁), n_mult*maximum(n₁))
            x_f, y_f, permutedims(n₁)
        end
        @series begin
            subplot := 2
            title := titles[2]
            color := cmapsim2
            seriestype := :heatmap
            n₂ = imag(sqrt.(ɛ-1im*sim.tls.D₀*F))
            clims := (-maximum(abs.(n₂)), maximum(abs.(n₂)))
            x_f, y_f, permutedims(n₂)
        end
        @series begin
            subplot := 3
            title := titles[3]
            color := cmapsim2
            clims := (-maximum(abs.(F)), maximum(abs.(F)))
            seriestype := :heatmap
            x_f, y_f, permutedims(F)
        end
        for j ∈ 1:3
            for k ∈ 1:2
                for i ∈ eachindex(sim.bnd.bc)
                    if isDirichlet(sim.bnd.bc[k][i])
                        linestyle := :solid
                    elseif isNeumann(sim.bnd.bc[k][i])
                        linestyle := :dash
                    end
                    if !isNone(sim.bnd.bl[k][i])
                        if isPMLout(sim.bnd.bl[k][i])
                            fill_color = :red
                        elseif isPMLin(sim.bnd.bl[k][i])
                            fill_color = :blue
                        end
                        @series begin
                            subplot := j
                            title := titles[j]
                            lw := 0
                            seriestype := :path
                            color := fill_color
                            alpha := .15
                            fill --> (0,.15,fill_color)
                            pml_data[i]
                        end
                    end
                    if isDirichlet(sim.bnd.bc[k][i]) | isNeumann(sim.bnd.bc[k][i])
                    @series begin
                        subplot := j
                        title := titles[j]
                        lw := 2
                        color := :black
                        seriestype := :path
                        boundary_data[i]
                    end
                end
            end
        end
    end

        for i ∈ 1:N
            ψ = ψ_plot[:,:,i]
            @series begin
                subplot := 3i+1
                clims --> (-maximum(abs.(ψ)), maximum(abs.(ψ)))
                color := cmapc
                x_t, y_t, permutedims(real(ψ))
            end
            @series begin
                subplot := 3i+2
                clims --> (-maximum(abs.(ψ)), maximum(abs.(ψ)))
                color := cmapc
                x_t, y_t, permutedims(imag(ψ))
            end
            @series begin
                subplot := 3i+3
                clims --> (0,maximum(abs2.(ψ)))
                color := cmapk
                x_t, y_t, permutedims(abs2.(ψ))
            end
        end
    else
        if iszero(N)
            layout --> (1,1)
            @series begin
                subplot := 1
                color := cmapsim1
                seriestype := :heatmap
                n = sqrt.(ɛ-1im*sim.tls.D₀*F)
                clims := (minimum(by.(n)), maximum(by.(n)))
                x_f, y_f, permutedims(by.(n))
            end
            for i ∈ eachindex(sim.bnd.bc)
                if isDirichlet(sim.bnd.bc[i])
                    linestyle := :solid
                elseif isNeumann(sim.bnd.bc[i])
                    linestyle := :dash
                end
                if !isNone(sim.bnd.bl[i])
                    if isPMLout(sim.bnd.bl[i])
                        fill_color = :red
                    elseif isPMLin(sim.bnd.bl[i])
                        fill_color = :blue
                    end
                    @series begin
                        subplot := 1
                        lw := 0
                        seriestype := :path
                        # color := fill_color
                        alpha := .15
                        fill --> (0,.15,fill_color)
                        pml_data[i]
                    end
                end
                if isDirichlet(sim.bnd.bc[i]) | isNeumann(sim.bnd.bc[i])
                    @series begin
                        subplot := 1
                        lw := 2
                        color := :black
                        seriestype := :path
                        boundary_data[i]
                    end
                end
            end
        else
            if by ∈ [abs, abs2]
                cmap = cmapk
                clims --> (0,1)
            else
                cmap = cmapc
                clims --> (-1,1)
            end
            if round(Int,N/3)==N/3
                N_col=3
                N_row = ceil(Int,N/3)
            elseif round(Int,N/2)==N/2
                N_col=2
                N_row = ceil(Int,N/2)
            elseif N==1
                N_col=1
                N_row=1
            else
                N_col = 3
                N_row = ceil(Int,N/3)
            end
            layout --> (N_row,N_col)
            for i ∈ 1:N
                ψ = ψ_plot[:,:,i]
                n₁ = real(sqrt.(ɛ-1im*sim.tls.D₀*F))
                renorm = (maximum(abs.(ψ)) - minimum(abs.(ψ)))/(maximum(n₁)-minimum(n₁))

                n₁ = n₁ .- minimum(n₁)
                n₁ = renorm*n₁
                @series begin
                    subplot := i
                    seriestype := :heatmap
                    seriesalpha := 1.
                    if by ∈ [abs, abs2]
                        color := cmapsim1
                    else
                        n₁ = n₁ .- maximum(abs.(ψ))
                        color := cmapsim1
                    end
                    x_f, y_f, permutedims(n₁)
                end
                @series begin
                    subplot := i
                    color := cmap
                    seriesalpha --> .85
                    # if by ∈ [abs, abs2]
                        # println(maximum(by.(ψ)))
                        # clims --> (0, maximum(by.(ψ)))
                    # else
                        # clims --> (-maximum(abs.(ψ)),+maximum(abs.(ψ)))
                    # end
                    x_t, y_t, permutedims(by.(ψ./maximum(abs.(ψ))))
                end
            end
        end
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


function prep_data(sim::Simulation{Tsys,Tbnd,Discretization{Cartesian}}, ψ, truncate::Bool) where Tsys<:System where Tbnd<:Boundary

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

    return (N, x_f_start, x_f_stop, x_f, x_t_start, x_t_stop, x_t[:],
        y_f_start, y_f_stop, y_f, y_t_start, y_t_stop, y_t[:],
        ε, F, ψ_plot,
        (boundary_data1_f, boundary_data2_f, boundary_data3_f, boundary_data4_f),
        (pml_data1, pml_data2, pml_data3, pml_data4) )
end


function prep_data(sim::Simulation{Tsys,Tbnd,Discretization{Polar}}, ψ, truncate::Bool) where Tsys<:System where Tbnd<:Boundary

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
    ε_itp = extrapolate(scale(interpolate(ε_temp,BSpline(Cubic(Periodic(OnGrid())))),r_f,θ_f),complex(NaN,NaN))
    ε = @. ε_itp(hypot(x_f,y_f'),sim.dis.x[2][1]+mod(atan(y_f',x_f)-sim.dis.x[2][1],2π))

    F_temp = Array{Float64}(undef, sim.dis.N[1]+2, sim.dis.N[2]+2)
    F_temp[2:end-1,2:end-1] = sim.sys.F
    @views F_temp[:,1] = F_temp[:,end-1]
    @views F_temp[:,end] = F_temp[:,2]
    @views F_temp[1,:] = F_temp[2,:]
    @views F_temp[end,:] = F_temp[end-1,:]
    F_itp = extrapolate(scale(interpolate(F_temp,BSpline(Cubic(Periodic(OnGrid())))),r_f,θ_f),NaN)
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

    return (N, x_f_start, x_f_stop, x_f, x_t_start, x_t_stop, x_t,
        y_f_start, y_f_stop, y_f, y_t_start, y_t_stop, y_t,
        ε, F, ψ_plot,
        (boundary_data1_f, boundary_data2_f, boundary_data3_f, boundary_data4_f),
        (pml_data1, pml_data2, pml_data3, pml_data4) )
end


################################################################################
########## ANIMATION
################################################################################
"""
    iterator = wave(sim, ψ; by=real, n=60, seriestype=:heatmap)

input for Plots.animate

Use cases:

`animate(wave(sim,ψ), file_name)` creates a .gif with filename

`animate(wave(sim,ψ; n=20), file_name, fps=10)` creates a 2 second movie

Note: default `fps`=20, and `n`=60, so default movie is 3 seconds long
"""
function wave(sim::Simulation, ψ; truncate=true, by=real, n=60)
    if 1 ∈ sim.dis.N
        return imap( ϕ->(sim, exp(-1im*ϕ)*ψ, by, truncate, 1), 0:2π/n:2π*(1-1/n))
    else
        return imap( ϕ->(sim, exp(-1im*ϕ)*ψ, by, truncate, 1, 1), 0:2π/n:2π*(1-1/n))
    end
end
