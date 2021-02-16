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

""" reset the model """
function reset!(model)
    fill!(model.c, 0)
    set_bc!(model)
end

function central_x(model)
    """ 
    Central difference scheme in
    x direction with periodic boundaries
    """
    return circshift(model.c, (0, -1)) + circshift(model.c, (0, 1))
end

function central_y(model)
    """ Central difference scheme
    in y direction
    """
    # inner_points = @view model.c[2:end-1, :]
    diff = circshift(model.c, (-1, 0)) + circshift(model.c, (1, 0))
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

function simulate(model, T; save=false, all=false)
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
    end
    
end

begin 
    """ simulate the model for T time units 

    if save=true and all=false: save only concentration along y-axis to file
    if save=true and all=true: save concentration at all grid points

    """
    # function simulate(model, T; save=false, all=false)
    #     if save
    #         if all
    #             return simulate_save_all(model, T)
    #         else
    #             return simulate_save(model, T)
    #         end
    #     else
    #         simulate_no_save(model, T)
    #     end
    # end

    """ Simulate for T time units and return concentration
    at all grid points to file """
    # function simulate_save_all(model, T)
    #     c_values = Array{Float64, 2}[]
    #     t = 0.0

    #     push!(c_values, model.c[:,:])#[:,1])
    #     while t < T
    #         step(model)
    #         push!(c_values, model.c[:,:])#[:,1])
    #         t += model.δt
    #     end
        
    #     return c_values
    # end

    # function simulate_save(model, T)
    #     c_values = Array{Float64}[]
    #     t = 0.0
    #     push!(c_values, model.c[:,1])
    #     while t < T
    #         step(model)
    #         push!(c_values, model.c[:,1])
    #         t += model.δt
    #     end
        
    #     return c_values
    # end

    # function simulate_no_save(model, T)
    #     t = 0.0
    #     while t < T
    #         step(model)
    #         t += model.δt
    #     end
    # end
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

    gif(anim, "anims/"*outpath*".gif", fps=10)
end

function plot_concentration_1D(filename, plotat; with_analytical=false)
    

    filename = "jlds/"*filename#*".jld"
    t_span = JLD.load(filename, "t")
    D = 1.0#JLD.load(filename, "D")
    plot_t_idxs = [argmin(abs.(t_span .- t)) for t in plotat]
    
    vals_to_plot = [JLD.load(filename, "$num")[:,1] for num in plot_t_idxs]
    t_vals = t_span[plot_t_idxs]
    println(t_vals)
    t_labels = reshape(["t=$(round(t, digits=4))" for t in t_vals], 1, :)
    yrange = range(0, stop=1, length=length(vals_to_plot[1]))

    colors = [:red, :blue, :yellow, :green, :purple]
    if with_analytical
        p = plot()
        for i ∈ 1:length(plotat)
            plot!(p, yrange, vals_to_plot[i], label=t_labels[i],
            xlabel="y", ylabel="Concentration", legend=:topleft)#, color=colors[i])

            analytical_y_range = range(0, stop=1, length=25)
            scatter!(p, analytical_y_range, c_analytical(analytical_y_range, plotat[i], D), 
                     label=false, color=:white)
        end
    else
        p = plot(yrange, vals_to_plot, label=t_labels,
        xlabel="y", ylabel="Concentration", legend=:topleft)
    end

    p
end


function plot_concentration_2D(filename, plotat)
    

    c_values, t_span = values(JLD.load("jlds/"*filename))
    N = size(c_values[1], 1)
    # t_span = range(0, stop=model.δt*length(c_values), length=length(c_values))
    plot_t_idx = argmin(abs.(t_span .- plotat))
    vals_to_plot = c_values[plot_t_idx]
    t = t_span[plot_t_idx]
    t_label = "t=$(round(t*1000, digits=2))"
    yrange = range(0, stop=1, length=N)

    heatmap(yrange, yrange, vals_to_plot, label=t_label,
                xlabel="x", ylabel="y")

end

erfc(arg) = (2/√π)*quadgk(t -> exp(-(t*t)), arg, Inf)[1]


""" Analytical solution to the diffusion equation """
function c_analytical(x, t, D, N=100)
    res = zeros(Float64, length(x))

    if iszero(t)
        """ If t=0.0, return zeros everywhere except at y=1"""
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
