using Plots

struct Concentration
    δt::Float64
    Γ::Float64
    c::Array{Float64, 2}
end

function initialize(N, D)

    δx = 1.0/N
    δy = 1.0/N
    δt = δx^2/(4.0*D)
    δt *= 0.99
    Γ = δt*D/δx^2


    c_array = zeros(Float64, N, N)
    model = Concentration(δt, Γ, c_array)
    set_bc!(model)
    return model

end

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

function finite_difference(model)
    """ Two-dimensional finite difference scheme """
end

function set_bc!(model)
    model.c[1,:] .= 0.0
    model.c[end,:] .= 1.0
end


function step(model)

    dcdx = central_x(model)
    dcdy = central_y(model)


    model.c .= model.c .+ model.Γ*(dcdx .+ dcdy - 4.0*model.c)
    set_bc!(model)

end

function animate(N, D, T, outpath; fps=10)

    model = initialize(N, D)
    N_t = round(Int, T/model.δt)

    anim = @animate for i ∈ 1:N_t
        t = i*model.δt
        heatmap(model.c, xlabel="x", ylabel="y", 
                title="t=$(round(t*1000, digits=2)) ms", 
                colorbar_title="Concentration")
        step(model)
    end

    gif(anim, outpath*".gif", fps=10)
end


function simulate(model, T; save=false)
    if save
        return simulate_save(model, T)
    else
        simulate_no_save(model, T)
    end
end

function simulate_save(model, T)
    c_values = Array{Float64}[]
    t = 0.0
    push!(c_values, model.c[:,1])
    while t < T
        step(model)
        push!(c_values, model.c[:,1])
        t += model.δt
    end
    
    return c_values
end

function simulate_no_save(model, T)
    t = 0.0
    while t < T
        step(model)
        t += model.δt
    end
end


function plot_concentration(N, D, T, plotat)
    

    model = initialize(N, D)
    c_values = simulate(model, T, save=true)
    t_span = range(0, stop=model.δt*length(c_values), length=length(c_values))
    plot_t_idxs = [argmin(abs.(t_span .- t)) for t in plotat]

    vals_to_plot = c_values[plot_t_idxs]
    t_vals = t_span[plot_t_idxs]
    t_labels = reshape(["t=$(round(t*1000, digits=2)) ms" for t in t_vals], 1, :)
    print(t_labels)
    plot(range(0, stop=1, length=N), vals_to_plot, label=t_labels,
    xlabel="y", ylabel="Concetration")


end

