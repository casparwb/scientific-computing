using Plots, JLD, QuadGK, PlotThemes
theme(:juno)

""" Concentration object """
struct Concentration
    δt::Float64            # time step size
    Γ::Float64             # stepping factor 
    c::Array{Float64, 2}   # concentration 
    t::Array{Float64}      # array of time values
end

""" Initialize a Concetration Object with N grid points and
diffusion coefficient D """
function initialize(N, D)

    δx = 1.0/N        # x-granularity
    δy = 1.0/N        # y-granularity
    δt = δx^2/(4.0*D) # (maximum) time step
    Γ = δt*D/δx^2     # stepping factor

    # initialize concentration to zero everywhere
    c_array = zeros(Float64, N, N)
    model = Concentration(δt, Γ, c_array, [0.0])
    set_bc!(model) # set boundary conditions
    return model

end

function set_bc!(model)
    model.c[1,:] .= 0.0
    model.c[end,:] .= 1.0
end


function step(model)

    dx = circshift(model.c, (0, -1)) + circshift(model.c, (0, 1))
    dy = circshift(model.c, (-1, 0)) + circshift(model.c, (1, 0))


    model.c[:,:] = model.c[:,:] .+ model.Γ*(dx .+ dy - 4.0*model.c[:,:])
    # model.t[1] += model.δt
    set_bc!(model)

end

""" Simulate the model for T time units """
function simulate(model, T)
    t = 0.0
    while t < T
        step(model)
        t += model.δt
    end
end

function save_to_file(N, D, T, outfile)
    

    model = initialize(N, D)
    # c_values = simulate(model, T, save=true, all=true)
    row = 2
    t = 0.0
    t_array = Float64[t]
    #sizehint!(t, round(Int, (T/model.δt)/saveevery))
    jldopen("jlds/"*outfile*".jld", "w") do file
        write(file, "1", model.c[:,:])
        while t < T
            step(model)
            write(file, "$(row)", model.c[:,:])
            row += 1
            t += model.δt
            push!(t_array, t)
        end
        write(file, "t", t_array)
        write(file, "D", D)
        write(file, "N", N)
    end
    
end



""" animate from file """
function animate_diffusion(filename, outpath; step=100, fps=10)

    filename = "jlds/"*filename*".jld" #  name of file to load
    file = jldopen(filename, "r")
    t_values = read(file, "t") #  array of time values 
    N_t = length(t_values)

    anim = @animate for i ∈ 1:step:N_t
        t = t_values[i]
        c = read(file, "$i")
        heatmap(c, xlabel="x", ylabel="y", 
                title="t=$(round(t, digits=2))", 
                colorbar_title="Concentration")
    end 
    close(file)
    gif(anim, "anims/"*outpath*".gif", fps=10)
end

function plot_concentration_1D(filename, plotat; outpath=nothing, with_analytical=false)
    
    # open file
    filename = "jlds/"*filename*".jld"
    file = jldopen(filename, "r")

    # get time-array and diffusion coefficient
    t_span = read(file, "t")
    D = read(file, "D")
    N = read(file, "N")

    # indices of given time spots to plot at
    plot_t_idxs = [argmin(abs.(t_span .- t)) for t in plotat]
    
    #vals_to_plot = [read(file, "$num")[:,1] for num in plot_t_idxs]
    t_vals = t_span[plot_t_idxs]
    t_labels = reshape(["t=$(round(t, digits=4))" for t in t_vals], 1, :)
    yrange = range(0, stop=1, length=N)

    colors = [:red, :blue, :yellow, :green, :purple]
    p = plot()
    for i ∈ 1:length(plotat)
        concentration = read(file, "$(plot_t_idxs[i])")[:,1]
        plot!(p, yrange, concentration, label=t_labels[i],
              xlabel="y", ylabel="Concentration", legend=:topleft)#, color=colors[i])
    end
    close(file)
    if with_analytical
            analytical_y_range = range(0, stop=1, length=25)
            for i ∈ 1:length(plotat)
                analytical_c = c_analytical(analytical_y_range, plotat[i], D)
            
                scatter!(p, analytical_y_range, analytical_c, 
                        label=false, color=:white)
            end
    end

    if !isnothing(outpath)
        savefig("plots/"*outpath*".svg")
    end
    p
end


function plot_concentration_2D(filename, plotat; outpath=nothing)
    

    # open file
    filename = "jlds/"*filename*".jld"
    file = jldopen(filename, "r")

    # get time-array and diffusion coefficient
    t_span = read(file, "t")
    D = read(file, "D")
    N = read(file, "N")


    # t_span = range(0, stop=model.δt*length(c_values), length=length(c_values))
    plot_t_idx = argmin(abs.(t_span .- plotat))[1]
    t = t_span[plot_t_idx]
    t_label = "t=$(round(t, digits=4))"
    yrange = range(0, stop=1, length=N)
    
    concentration = read(file, "$plot_t_idx")
    heatmap(yrange, yrange, concentration, label=t_label,
            xlabel="x", ylabel="y", cbar_title="Concentration",
            title=t_label)

    if !isnothing(outpath)
        savefig("plots/"*outpath*".svg")
    end
end

erfc(arg) = (2/√π)*quadgk(t -> exp(-(t*t)), arg, Inf)[1]


""" Analytical solution to the diffusion equation """
function c_analytical(x, t, D, N=100)

    res = zeros(Float64, length(x)) # array for storing result

    """ If t=0.0, return zeros everywhere except at y=1"""
    if iszero(t)
        res[end] = 1.0
        return res
    end    

    denom = 1.0/(2*√(D*t)) # denominator factor
    for k ∈ 1:length(x)    # x-values
        for i ∈ 0:N        # sum for each x-value
            arg1 = (1 - x[k] + 2*i)*denom
            arg2 = (1 + x[k] + 2*i)*denom 

            res[k] += (erfc(arg1) - erfc(arg2))
        end
    end

    res
end
