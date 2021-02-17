using Plots, PlotThemes, QuadGK
theme(:juno)


struct Concentration
    c::Array{Float64, 2}
    ω::Float64
    idxs::Union{Dict{Tuple{Int8, Int8}, Int8}}
end

""" Set boundary conditions in y-direction (source and sink)"""
function set_bc!(array)
    array[1,:] .= 0.0
    array[end,:] .= 1.0
end

function initialize(N, ω=1.0; objects=nothing)
    """ Objects: array of tuples with (x₀, x₁, y₀, y₁)"""
    δx = 1.0/N
    δy = 1.0/N

    c_array = zeros(Float64, N, N)
    set_bc!(c_array)

    all_idxs = 0:N+1
    idx_dict = Dict{Tuple{Int8, Int8}, Int8}([(Int8(y), Int8(x)) => 1 for x in all_idxs, y in all_idxs])

    if isnothing(objects)
        model = Concentration(c_array, ω, idx_dict)#, c_array[:,:])
        return model
    else
        axis_range = range(0, stop=1, length=N)
        for obj in objects
            x₀, x₁, y₀, y₁ = obj
            x_idxs = findall(x₀ .<= axis_range .<= x₁)
            y_idxs = findall(y₀ .<= axis_range .<= y₁)

            for x in x_idxs, y in y_idxs
                idx_dict[(Int8(y), Int8(x))] = 0
            end
            
        end

        model = Concentration(c_array, ω, idx_dict)
        return model
    end

end

function step_Jacobi!(model, ϵ)
    flag = true
    c_old = model.c[:,:]
    dx = circshift(model.c, (0, -1)) + circshift(model.c, (0, 1))
    dy = circshift(model.c, (-1, 0)) + circshift(model.c, (1, 0))

    @. model.c[:,:] = 0.25*(dx + dy)
    set_bc!(model.c)
    if maximum(abs.(model.c .- c_old)) < ϵ 
        flag = false
    end

    return flag

end


function step_GaussSeidel!(model, ϵ)
    
    ny, nx = size(model.c)

    flag = true
    maxdiff = 0
    for i = 2:ny-1 
        old_row = model.c[i,:]
        for j = 1:nx
            up = model.c[i+1, j]
            down = model.c[i-1, j]

            if j == 1
                left = model.c[i, end]
                right = model.c[i, j+1]
            elseif j == nx
                left = model.c[i, end-1]
                right = model.c[i, 1]
            else
                left = model.c[i, j-1]
                right = model.c[i, j+1]
            end                

            model.c[i, j] = 0.25*(up + 
                                 down +
                                 right +
                                 left)
        end
        rowdiff = maximum(abs.(old_row .- model.c[i,:]))
        if 0 < rowdiff > maxdiff
            maxdiff = rowdiff
        end
    end
    if maxdiff < ϵ
        flag = false
    end

    set_bc!(model.c)
    return flag
end

function step_SOR!(model, ϵ)
    ny, nx = size(model.c)

    flag = true
    maxdiff = 0

    for i = 2:ny-1
        old_row = model.c[i,:]
        for j = 1:nx
            iszero(model.idxs[Int8.((i, j))]) && continue

            up = model.c[i+1, j]
            down = model.c[i-1, j]
            center = model.c[i, j]

            if j == 1
                left = model.c[i, end]
                right = model.c[i, j+1]
            elseif j == nx
                left = model.c[i, end-1]
                right = model.c[i, 1]
            else
                left = model.c[i, j-1]
                right = model.c[i, j+1]
            end                

            model.c[i, j] = 0.25*model.ω*(up + 
                                 down +
                                 right +
                                 left) + (1 - model.ω)*center
        end
        rowdiff = maximum(abs.(old_row .- model.c[i,:]))
        if 0 < rowdiff > maxdiff
            maxdiff = rowdiff
        end
    end
    if maxdiff < ϵ
        flag = false
    end

    set_bc!(model.c)
    return flag
end

function step_SOR_insulating!(model, ϵ)
    ny, nx = size(model.c)

    flag = true
    maxdiff = 0

    for i = 2:ny-1
        old_row = model.c[i,:]
        for j = 1:nx
            iszero(model.idxs[Int8.((i, j))]) && continue # inside the object: c = 0
            if iszero(model.idxs[Int8.((i+1, j))]) 
                center = model.c[i, j]
                up = center
                down = model.c[i-1, j]
                left = model.c[i, j-1]
                right = model.c[i, j+1]
            elseif iszero(model.idxs[Int8.((i-1, j))])
                center = model.c[i, j]
                up = model.c[i+1, j]
                down = center
                left = model.c[i, j-1]
                right = model.c[i, j+1]
            elseif iszero(model.idxs[Int8.((i, j+1))])
                center = model.c[i, j]
                up = model.c[i+1, j]
                down = model.c[i-1, j]
                left = model.c[i, j-1]
                right = center
            elseif iszero(model.idxs[Int8.((i, j-1))])
                center = model.c[i, j]
                up = model.c[i+1, j]
                down = model.c[i-1, j]
                left = center
                right = model.c[i, j+1]
            else
                up = model.c[i+1, j]
                down = model.c[i-1, j]
                center = model.c[i, j]

                if j == 1
                    left = model.c[i, end]
                    right = model.c[i, j+1]
                elseif j == nx
                    left = model.c[i, end-1]
                    right = model.c[i, 1]
                else
                    left = model.c[i, j-1]
                    right = model.c[i, j+1]
                end                
            end

            model.c[i, j] = 0.25*model.ω*(up + 
                                 down +
                                 right +
                                 left) + (1 - model.ω)*center
        end
        rowdiff = maximum(abs.(old_row .- model.c[i,:]))
        if 0 < rowdiff > maxdiff
            maxdiff = rowdiff
        end
    end
    if maxdiff < ϵ
        flag = false
    end

    set_bc!(model.c)
    return flag
end



function simulate(model, solver, ϵ=1e-8)


    if solver == "Jacobi"
        stepfunc = step_Jacobi!
    elseif solver =="GaussSeidel"
        stepfunc = step_GaussSeidel!
    elseif solver == "SOR_Insulating"
        stepfunc = step_SOR_insulating!
    else
        stepfunc = step_SOR!
    end


    counter = 0
    flag = true
    while flag
        flag = stepfunc(model, ϵ)
        counter += 1
    end

    return counter
end



""" Analytical solution to the diffusion equation """
function c_analytical(x, t, D, N=100)
    res = zeros(Float64, length(x))
    erfc(arg) = (2/√π)*quadgk(t -> exp(-(t*t)), arg, Inf)[1]

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

function compare_solvers(ω, N=50, ϵ=1e-6)

    model_Jacobi = initialize(N)
    model_GaussSeidel = initialize(N)
    model_SOR = initialize(N, ω)
    
    k_Jacobi = simulate(model_Jacobi, "Jacobi", ϵ)
    k_GaussSeidel = simulate(model_GaussSeidel, "GaussSeidel", ϵ)
    k_SOR = simulate(model_SOR, "SOR", ϵ)

    y_range = range(0, stop=1.0, length=N)

    y_Jacobi = model_Jacobi.c[:,1]
    y_GaussSeidel = model_GaussSeidel.c[:,1]
    y_SOR = model_SOR.c[:,1]

    y_analytical = c_analytical(y_range, 1.0, 1.0, 50)


    p1 = plot(y_range, y_Jacobi,      label="Jacobi", xlabel="y", ylabel="Concentration")
    plot!(p1, y_range, y_GaussSeidel, label="Gauss-Seidel")
    plot!(p1, y_range, y_SOR,         label="SOR")
    scatter!(p1, y_range, y_analytical,  label="Analytical", markersize=3.5 )

    δ_values = 10 .^ range(-2, stop=-8, length=20)
    ω_values = [1.75, 1.8, 1.95]
    results = zeros(length(δ_values), 5)

    solvers = ["Jacobi" "GaussSeidel" "SOR (ω=1.75)" "SOR (ω=1.8)" "SOR (ω=1.95)"]

    for i = 1:length(δ_values)
        model_Jacobi = initialize(N)
        model_GaussSeidel = initialize(N)
        model_SOR1 = initialize(N, ω_values[1])
        model_SOR2 = initialize(N, ω_values[2])
        model_SOR3 = initialize(N, ω_values[3])

        k_Jacobi = simulate(model_Jacobi, "Jacobi", δ_values[i])
        k_GaussSeidel = simulate(model_GaussSeidel, "GaussSeidel", δ_values[i])
        k_SOR1 = simulate(model_SOR1, "SOR", δ_values[i])
        k_SOR2 = simulate(model_SOR2, "SOR", δ_values[i])
        k_SOR3 = simulate(model_SOR3, "SOR", δ_values[i])

        results[i, 1] = k_Jacobi
        results[i, 2] = k_GaussSeidel
        results[i, 3] = k_SOR1
        results[i, 4] = k_SOR2
        results[i, 5] = k_SOR3
    end

    p2 = scatter(δ_values, results, xscale=:log10, yscale=:log10,
                 label=solvers, xlabel="log₁₀(δ)", ylabel="log₁₀(# iterations)")
    

    savefig(p1, "model_analytical.svg")
    savefig(p2, "delta_vs_iterations.svg")
end

function find_optimal_ω(ϵ=1e-8; objects=nothing)

    ω_values = range(1.71, stop=1.99, length=100)
    N_values = range(10, stop=100, length=10)
    ks = zeros(Int, length(N_values), length(ω_values))

    for i = 1:length(N_values)
        for j = 1:length(ω_values)
            model_SOR = initialize(Int(N_values[i]), ω_values[j], objects=objects)
            ks[i, j] = simulate(model_SOR, "SOR", ϵ)
            println("Finished $((i, j))")
        end
    end


    best_ω = ω_values[argmin(ks)]
    if isnothing(objects)
        println("Optimal ω-values w/o objects: $(best_ω)")
    else
        println("Optimal ω-values with objects: $(best_ω)")
    end
        

    return ω_values, N_values, ks
end

function plot_optimal_w(filename)

    w, N, k, o = JLD.load(filename*".jld", "omega", "N", "k", "objects")
    if isnothing(o)  o = [] end

    heatmap(w, N, log10.(k), xlabel="ω", ylabel="N", cbar_title="log₁₀(# iterations)")
    title!("# sink objects: $(length(o))")
    savefig("$filename"*"_heatmap.svg")

    plot(w, k[5, :], xlabel="ω", ylabel="# iterations", label="N=$(N[5])", yscale=:log10)
    title!("# sink objects: $(length(o))")
    savefig("$filename"*"_n=50.svg")

    plot(w, N, log10.(k), st=:surf, xlabel="ω", ylabel="N", zlabel="log₁₀(# iterations)")
    title!("# sink objects: $(length(o))")
    savefig("$filename"*"_surf.svg")

end

function animate_ti_diffusion(model, outpath; solver="SOR", n=500, step=10, fps=10, ϵ=1e-8)

    axis = range(0, stop=1, length=size(model.c, 1))
    #heatmap(axis, axis, model.c, xlabel="x", ylabel="y", 
    #            colorbar_title="Concentration", title="i = 0")

    if solver == "Jacobi"
        stepfunc = step_Jacobi!
    elseif solver =="GaussSeidel"
        stepfunc = step_GaussSeidel!
    elseif solver == "SOR_Insulating"
        stepfunc = step_SOR_insulating!
    else
        stepfunc = step_SOR!
    end

    anim = @animate for i ∈ 1:n
        flag = stepfunc(model, ϵ)
        heatmap(axis, axis, model.c, xlabel="x", ylabel="y", 
                colorbar_title="Concentration", title="k = $i")
        println("$(round(i/n*100, digits=2)) %")
    end when flag

    gif(anim, "anims/"*outpath*".gif", fps=fps)
end


function plot_concentration_2D(model, ϵ = 1e-8; at_k = nothing, title=nothing, solver=step_SOR!)

    if isnothing(at_k)
        p = heatmap(model.c, xlabel="x", ylabel="y", cbar_title="Concentration",
                    title=title)
    else
        counter = 0
        for i ∈ 1:at_k
            solver(model, ϵ)
        end

        p = heatmap(model.c, xlabel="x", ylabel="y", cbar_title="Concentration",
                    title="k = $at_k")

    end

    return p
end