using Plots, PlotThemes, StatsBase
theme(:juno)

""" Concentration object """
struct Concentration
    c::Array{Float64, 2} # matrix of concentration in domain
    ω::Float64           # relaxation factor 
    η::Float64           # 
    N::Int64
    idxs::Union{Dict{Tuple{Int8, Int8}, Int8}} # dictionary of indices (for identifying where objects are)
end

struct Concentration2
    c::Array{Float64, 2} # matrix of concentration in domain
    ω::Float64           # relaxation factor 
    η::Float64           # 
    N::Int64
    idxs::Union{Dict{Tuple{Int8, Int8}, Int8}} # dictionary of indices (for identifying where objects are)
end

""" Set boundary conditions in y-direction (source and sink)"""
function set_bc!(array)
    array[1,:] .= 0.0
    array[end,:] .= 1.0
end

""" Initialize a concentration object,
    starting with a linear concentration gradient
    - objects: array of tuples with (x₀, x₁, y₀, y₁) where objects should be
"""
# function initialize(N, ω=1.95; objects=nothing)
#     δx = 1.0/N
#     δy = 1.0/N

#     c_array = zeros(Float64, N, N) # initialize array 
#     for i in range(N)
#         c_array[:,i] .= range(0, stop=1, length=N)
#     end

#     ### Set up index matrix ###
#     all_idxs = 0:N+1 

#     # set all indices to 1 initially
#     idx_dict = Dict{Tuple{Int8, Int8}, Int8}([(Int8(y), Int8(x)) => 1 for x in all_idxs, y in all_idxs])

#     if isnothing(objects) # if no objects, just return object
#         model = Concentration(c_array, ω, idx_dict)
#         return model
#     else # if object, set indices where object lies to 0
#         axis_range = range(0, stop=1, length=N)
#         for obj in objects
#             x₀, x₁, y₀, y₁ = obj
#             x_idxs = findall(x₀ .<= axis_range .<= x₁)
#             y_idxs = findall(y₀ .<= axis_range .<= y₁)

#             for x in x_idxs, y in y_idxs
#                 idx_dict[(Int8(y), Int8(x))] = 0 # set indices to 0
#             end
            
#         end
        
#         model = Concentration(c_array, ω, idx_dict)
#         return model
#     end

# end

function initialize(N, η, ω=1.95; seed=N/2)
    δx = 1.0/N
    δy = 1.0/N

    c_array = zeros(Float64, N, N) # initialize array 
    for i =1:N
        c_array[:,i] .= range(0, stop=1, length=N)
    end

    all_idxs = 1:N 

    # set all indices to 1 initially
    idx_dict = Dict{Tuple{Int8, Int8}, Int8}([(Int8(y), Int8(x)) => 1 for x in all_idxs, y in all_idxs])
    seed_idx = (Int8(1), Int8(N/2))
    idx_dict[seed_idx] = 0 # growth are characterized by 0

    # candidates are characterized by 2
    initial_candidates = [(1, seed_idx[2]+1), (1, seed_idx[2]-1)]
    idx_dict[initial_candidates[1]] = 2
    idx_dict[initial_candidates[2]] = 2
        
    c_array[seed_idx...] = 0
    model = Concentration2(c_array, ω, η, N, idx_dict)
    return model


end

function choose_candidate(model)
    candidate_idxs = findall(x -> x == 2, model.idxs) # current possible growth candidates

    concentrations = model.c[ CartesianIndex.(candidate_idxs)] .^ model.η # c^η   
    probs = concentrations / sum(concentrations) # c^η/sum(c^η)

    choice = sample(candidate_idxs, ProbabilityWeights(probs)) # new growth cell
    neighbours = [(choice[1]+1, choice[2]),
                  (choice[1]-1, choice[2]),
                  (choice[1], choice[2]+1),
                  (choice[1], choice[2]-1)] # neigbours of new cell

    # only choose neighbours who are not outside the domain
    new_candidate_candidates = [idx for idx in neighbours if !(any(any.((idx .== 0, idx .== model.N+1))))]

    for new_candidate_candidate in new_candidate_candidates
        if model.idxs[new_candidate_candidate] == 1 # if not already parth of growth, set as new candidate
             model.idxs[new_candidate_candidate] = 2
        end
    end

    model.idxs[choice] = 0 # set index identity to 0
    return choice
end

function step!(model, ϵ=1e-6; diffuse=false)

    growth_choice = choose_candidate(model) # growth choice
    model.c[CartesianIndex(growth_choice)] = 0 # set concentration here to 0
    if diffuse
        diffuse!(model, ϵ)
    end
end

function simulate(model, N; diffuse_every=1)

    for i = 1:N
        diffuse = false
        if iszero(i%diffuse_every)
            diffuse = true
        end

        step!(model, diffuse=diffuse)
    end

end

""" Diffuse concentration fields one time """
function diffuse!(model, ϵ=1e-6)
    ny, nx = size(model.c)

    flag = true
    maxdiff = 0

    for i = 2:ny-1
        old_row = model.c[i,:]
        for j = 1:nx
            # if the current indices are inside the object domain, just skip the iteration (leave them to be zero)
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

function plot_growth(model)

    domain = ones(Int8, size(model.c))

    growth_idxs = CartesianIndex.(findall(x -> x == 0, model.idxs))
    candidate_idxs = CartesianIndex.(findall(x -> x == 2, model.idxs))

    domain[growth_idxs] .= 0
    domain[candidate_idxs] .= 2

    ax = range(0, stop=1, length=model.N)

    heatmap(ax, ax, domain, xlabel="x", ylabel="y", title="Growth", colorbar=false)
end

function plot_concentration(model)


    ax = range(0, stop=1, length=model.N)

    heatmap(ax, ax, model.c, xlabel="x", ylabel="y", colorbar_title="Concentration",
            title="Concentration Field", clims=(0, 1))
end

function plot_both(model)

    p1 = plot_growth(model)
    p2 = plot_concentration(model)

    plot(p1, p2, size=(1000, 400))

end

function animate_growth(model, N; step=1)

    anim = @animate for i = 1:N
        p1 = plot_growth(model)
        p2 = plot_concentration(model)
        plot(p1, p2, size=(1000, 400))
        step!(model, 1e-6, diffuse=true)
    end every step

    gif(anim, "test2.gif")
end