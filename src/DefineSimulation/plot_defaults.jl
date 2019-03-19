function __init__()
    if !haskey(ENV, "SCALAR_FDFD_COLOR_THEME")
        SCALAR_FDFD_COLOR_THEME = :default
    else
        SCALAR_FDFD_COLOR_THEME = Symbol(ENV["SCALAR_FDFD_COLOR_THEME"])
    end

    global COLOR_SCHEME = SCALAR_FDFD_COLOR_THEME
    set_colors_by_scheme()
end


function set_colors_by_scheme()
    if COLOR_SCHEME == :default
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :dark
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :ggplot2
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :juno
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :lime
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :orange
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :sand
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :solarized
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :solarized_light
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :wong
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    elseif COLOR_SCHEME == :wong2
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    else
        global BAND_COLOR = :lightgrey
        global BAND_WIDTH = 1.2
        global BAND_STYLE = :solid
        global GAP_COLOR = :lightgrey
        global GAP_WIDTH = 1
        global GAP_STYLE = :dash
        global DISPERSION_WIDTH = 3
        global DISPERSION_STYLE = :solid
    end
    return nothing
end


"""
    cmapc, cmapk, cmapsim1, cmapsim2, n_mult, F_sign = fix_colormap(theme)
"""
function fix_colormap(theme)

    F_sign = +1
    n_mult=1

    if theme ∈ [:dark, :juno]
        cmapc=:bkr
        cmapk=:ice
        cmapsim1=:dimgray
        cmapsim2=:bky
        n_mult=1.3
    elseif theme ∈ [:solarized]
        cmapc=:bky
        cmapk=:solar
        cmapsim1=:viridis
        cmapsim2=:bky
        n_mult=1.3
    elseif theme ∈ [:orange]
        cmapc=:bky
        cmapk=:haline
        cmapsim1=:inferno
        cmapsim2=:bky
        n_mult=1.1
    elseif theme ∈ [:lime]
        cmapc=:bky
        cmapk=:haline
        cmapsim1=:inferno
        cmapsim2=:bky
        n_mult=1.0
    elseif theme ∈ [:solarized_light]
        cmapc=:YlOrBr
        cmapk=:Greys
        cmapsim1=:YlOrBr
        cmapsim2=:RdYlBu
        F_sign = -1
    else
        cmapc=:RdBu
        cmapk=:Greys
        cmapsim1=:Greys
        cmapsim2=:RdGy
        F_sign = -1
    end

    return cmapc, cmapk, cmapsim1, cmapsim2, n_mult, F_sign
end
