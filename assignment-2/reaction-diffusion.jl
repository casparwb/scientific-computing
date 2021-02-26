using Plots, PlotThemes
theme(:juno)


struct ReactionDiffusionModel
    Dᵤ::Real
    Dᵥ::Real
    f::Real
    k::Real
    δx::Real
    δt::Real
    u::Array{Float64, 2}
    v::Array{Float64, 2}
end

function initialize(N; δt=1, δx=1, Dᵤ=0.16, Dᵥ=0.08, f=0.35, k=0.06, u₀=0.5, v₀=0.25)

    u = ones(N, N)*u₀
    v = zeros(N, N)

    Δ = round(Int, N/50)

    v[Int(N/2)-Δ:Int(N/2)+Δ, Int(N/2)-Δ:Int(N/2)+Δ] .= v₀

    return ReactionDiffusionModel(Dᵤ, Dᵥ, f, k, δx, δt, u, v)
end

""" Step model once using Successive Over-Relaxation scheme 
with periodic boundary conditions in both dimensions"""
function ∇²(array, ϵ=1e-8, ω=1.92)
    ny, nx = size(array)

    flag = true
    maxdiff = 0

    result = similar(array)
    for i = 1:ny
        old_row = array[i,:]
        for j = 1:nx
            # if the current indices are inside the object domain, just skip the iteration (leave them to be zero)
            #iszero(model.idxs[Int8.((i, j))]) && continue 

            center = array[i, j]

            if i == 1
                up = array[i+1, j]
                down = array[ny, j]
            elseif i == ny
                up = array[1, j]
                down = array[end-1, j]
            else
                up = array[i+1, j]
                down = array[i-1, j]
            end        

            if j == 1
                left = array[i, end]
                right = array[i, j+1]
            elseif j == nx
                left = array[i, end-1]
                right = array[i, 1]
            else
                left = array[i, j-1]
                right = array[i, j+1]
            end                

            result[i, j] = 0.25*ω*(up + 
                                 down +
                                 right +
                                 left) + (1 - ω)*center
        end
        rowdiff = maximum(abs.(old_row .- array[i,:]))
        if 0 < rowdiff > maxdiff
            maxdiff = rowdiff
        end
    end
    if maxdiff < ϵ
        flag = false
    end

    #set_bc!(array)
    return result
end

function ∇²_2(array)

    result = similar(array)

    dx = circshift(array, (0, -1)) + circshift(array, (0, 1)) # cᵏᵢ,ⱼ₊₁ + cᵏᵢ,ⱼ₋₁
    dy = circshift(array, (-1, 0)) + circshift(array, (1, 0)) # cᵏᵢ₊₁,ⱼ + cᵏᵢ₋₁,ⱼ

    result[:,:] = (dx .+ dy .- 4.0*array[:,:]) # step the model

    return result

end

function gray_scott!(model)

    ∇²u = ∇²_2(model.u)/model.δx^2
    ∇²v = ∇²_2(model.v)/model.δx^2

    uv² = @. model.u * model.v * model.v

    ∂u∂t = @. model.Dᵤ*∇²u - uv² + model.f*(1 - model.u)
    ∂v∂t = @. model.Dᵥ*∇²v - uv² - (model.f + model.k)*model.v

    model.u[:,:] = model.u[:,:] .+ model.δt*∂u∂t
    model.v[:,:] = model.v[:,:] .+ model.δt*∂v∂t

end



function plot_concentrations(model)

    # cmax = max(maximum(model.u), maximum(model.v))
    # cmin = min(minimum(model.u), minimum(model.v))

    p1 = heatmap(model.u, title="U", ylabel="y")#, clims=(0, 1))
    p2 = heatmap(model.v, title="V", xlabel="x", ylabel="y", colorbar_title="Concentration")
    plot(p1, p2, layout=(2, 1), size=(600, 800))

end

function animate_concentrations(n_steps, model)


    anim = @animate for i = 1:n_steps
        gray_scott!(model)
        plot_concentrations(model)
    end

    gif(anim, "test.gif")
end

