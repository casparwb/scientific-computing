using Plots, PlotThemes
theme(:juno)


struct ReactionDiffusionModel
    Dᵤ::Real       # U-concentration diffusion coefficient
    Dᵥ::Real       # V-concentration diffusion coefficient
    f::Real        # U supply rate
    k::Real        # V decay rate
    δx::Real       # Grid granulatity
    δt::Real       # Temporal granularity
    u::Array{Float64, 2} # U-concentration array
    v::Array{Float64, 2} # V-concentration array
end

function initialize(N; δt=1, δx=1, Dᵤ=0.16, Dᵥ=0.08, f=0.035, k=0.06, u₀=0.5, v₀=0.25,ξ=0.1)

    u = ones(N, N)*u₀
    v = zeros(N, N)

    Δ = round(Int, N/50) # size of square where initial concentration is v₀
    v[Int(N/2)-Δ:Int(N/2)+Δ, Int(N/2)-Δ:Int(N/2)+Δ] .= v₀
    v[:,:] .+= ξ*rand(N, N) # add noise

    return ReactionDiffusionModel(Dᵤ, Dᵥ, f, k, δx, δt, u, v)
end

""" Compute Laplacian of array with either zero- or periodic boundaries """
function ∇²(array; boundaries="periodic")

    ny, nx = size(array)
    result = zeros(ny, nx)

    if boundaries == "periodic"
        dx = circshift(array, (0, -1)) + circshift(array, (0, 1)) # aᵏᵢ,ⱼ₊₁ + aᵏᵢ,ⱼ₋₁
        dy = circshift(array, (-1, 0)) + circshift(array, (1, 0)) # aᵏᵢ₊₁,ⱼ + aᵏᵢ₋₁,ⱼ
        result[:,:] = @. (dx + dy - 4.0*array[:,:]) # step the model
    else
        for i = 2:ny-1
            for j = 2:nx-1
                result[i, j] = array[i, j-1] + 
                array[i, j+1] +
                array[i-1, j] + 
                array[i+1, j] -
                4.0*array[i, j]
            end
        end
    end
    
    return result

end

""" Step the Gray-Scott model one time step """
function gray_scott!(model)

    # ∇²u = ∇²(model.u, boundaries="zero")/model.δx^2
    # ∇²v = ∇²(model.v, boundaries="zero")/model.δx^2

    ∇²u = ∇²(model.u)/model.δx^2
    ∇²v = ∇²(model.v)/model.δx^2

    uv² = @. model.u * model.v * model.v

    ∂u∂t = @. model.Dᵤ*∇²u - uv² + model.f*(1.0 - model.u)
    ∂v∂t = @. model.Dᵥ*∇²v + uv² - (model.f + model.k)*model.v

    model.u[:,:] .= model.u[:,:] .+ model.δt*∂u∂t
    model.v[:,:] .= model.v[:,:] .+ model.δt*∂v∂t

end


""" Plot both U and V concentrations """
function plot_concentrations(model; outfile=nothing)

    ax = range(0, stop=1, length=size(model.u, 1))

    p1 = heatmap(ax, ax, model.u, title="U", ylabel="y", xlabel="x", cbar=false)#)
    p2 = heatmap(ax, ax, model.v, title="V", xlabel="x", cbar=false)#, colorbar_title="Concentration")
    plot(p1, p2, layout=(1, 2), size=(800, 400))

    if !(isnothing(outfile))
        savefig("figures/"*outfile*".svg")
    end
end

""" Plot only one concentration """
function plot_concentration(model, c="V"; outfile=nothing, title=nothing)
    if c == "V"
        array = model.v
    else
        array = model.v
    end

    ax = range(0, stop=1, length=size(model.u, 1))

    p = heatmap(ax, ax, array, xlabel="x", ylabel="y", cbar=false, aspect_ratio=1)
    !(isnothing(title)) && title!(title)
    if !(isnothing(outfile))
        savefig("figures/"*outfile*".svg")   
    end 

    return p
end


function animate_concentrations(n_steps, model, outfile; step=10)

    anim = @animate for i = 1:n_steps
        gray_scott!(model)
        plot_concentrations(model)
    end every step

    gif(anim, "anims/"*outfile*".gif")
end

function simulate!(model, T)
    
    t = 0.0
    while t < T
        gray_scott!(model)
        t += model.δt
    end

end


function different_parameters(;N=100, T=5000)

    f_values = [0.02, 0.035, 0.05]
    k_values = [0.03, 0.05, 0.06]

    for f in f_values
        for k in k_values
            model = initialize(N, f=f, k=k)
            simulate!(model, T)
            plot_concentration(model, outfile="f=$(f)_k=$k", title="f = $f, k = $k")
        end
    end

end

# animate_concentrations(5000, initialize(256), "f=0.035_k=0.06")