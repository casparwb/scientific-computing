### Set 1: Vibrating String ###
using Plots, PlotThemes
theme(:juno)

""" Object for modeling the string """
struct VibratingString
    c::Float64   # wave velocity
    Δx::Float64  # x-grid size
    Δt::Float64  # t-grid size
    Γ::Float64   # step parameter
    t::Array{Float64} # array of time values
    Ψ::Array{Float64, 2} # array of amplitude values
end

""" set boundary conditions for ψ(tᵢ) and Ψ(tᵢ₊₁) """
function set_bc!(string)
    string.Ψ[[1, end], :] .= 0.0
end

""" 
Initialize and return String object

Input
------
    - init_func: function(::array), function for initializing Ψ(x) at t=0
    - N: Int, number of grid points
    - Δt: Number, time step size
    
Returns
------
    - String object 
"""
function initialize(init_func, N, Δt)

    L = 1.0
    c = 1.0

    Δx = L/N
    x = range(0, stop=L, length=N)

    Ψ = zeros(Float64, N, 2)
    Ψ[:,1] = init_func(x)
    Ψ[:,2] = Ψ[:,1]

    Γ = (Δt*c/Δx)^2
    string = VibratingString(c, Δx, Δt, Γ, [0.0], Ψ)

    set_bc!(string)
    return string
end

""" Step String model one time step. In-place operation."""
function step!(model)

    term1 = circshift(model.Ψ[:,2], -1) - # Ψᵢ₊₁,ⱼ
            2*model.Ψ[:,2]              + # Ψᵢ,ⱼ
            circshift(model.Ψ[:,2], 1)    # Ψᵢ₋₁,ⱼ

    term2 = 2*model.Ψ[:,2] - model.Ψ[:,1] # Ψᵢ,ⱼ - Ψᵢ,ⱼ₋₁
    model.Ψ[:,1] = model.Ψ[:,2] # swap new state with old state
    model.Ψ[:,2] = model.Γ*term1 + term2  # update Ψᵢ,ⱼ₊₁ 
    set_bc!(model) # set boundary conditions

end


""" Simulate for T time units """
function simulate(model, T)
    t = 0.0
    while t < T
        step!(model)
        t += model.Δt
        push!(model.t, t)
    end
end


function animate_string(model, T, outpath; step=1, fps=1)

    N = Int(T/model.Δt)
    t = 0.0
    anim = @animate for i ∈ 1:N
        step!(model)
        t += model.Δt
        plot(model.Ψ[:,2], ylims=(-1, 1), title="t=$(round(t, digits=2))", 
                legend=false, xlabel="x", ylabel="Ψ")
    end every step

    gif(anim, outpath*".gif", fps = fps)
end


function plot_wave(model; title=nothing, label=nothing)

    x = range(0, stop=1, length=size(model.Ψ, 1))

    p = plot(x, model.Ψ[:,2], label=label)
    xlabel!("x")
    ylabel!("Ψ")
    ylims!(-1, 1)
    if !isnothing(title)
        title!(title)
    end
    p
end

function visualize(N, Δt, T; fps=10, nplots=10)

    f1(x) = @. sin(2π*x)
    f2(x) = @. sin(5π*x)
    function f3(x) 
        out = zeros(length(x))
        where_inner = findall(0.2 .< x .< 0.4)
        out[where_inner] = f2(x[where_inner])
        return out
    end

    init_funcs = [f1, f2, f3]
    outnames = ["2pi", "5pi", "inner_5pi"]

    x = range(0, stop=1, length=N)

    t_plots = range(0, stop=T, length=nplots) # time at which to plot
    for (init, outname) in zip(init_funcs, outnames)
        ### plot ###
        p = plot()
        string = initialize(init, N, Δt)
        plot!(p, x, string.Ψ[:,2], title=outname, label="t = $(round(t_plots[1], digits=2))")
        for t_plot in t_plots[2:end]
            simulate(string, t_plot)
            plot!(p, x, string.Ψ[:,2], label="t = $(round(t_plot, digits=2))")
        end
        xlabel!("x")
        ylabel!("Ψ")
        savefig("plots/"*outname*".svg")


        ### animate ###
        animate_string(initialize(init, N, Δt), T, "anims/"*outname, step=10, fps=fps)
    end


end
