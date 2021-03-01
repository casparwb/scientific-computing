include("code.jl")
using Plots, PlotThemes#, StatsBase
theme(:juno)


const steps = ["up", "down", "left", "right"]
const moves = Dict("up" => [1, 0],
                   "down" => [-1, 0],
                   "left" => [0, -1],
                   "right" => [0, 1])

struct RandomWalker
    pos::Array{Int8, 1}
end

struct Growth
    indices::Union{Dict{Tuple{Int8, Int8}, Int8}}
end

struct Cluster
    N::Int16
    domain::Array{Int8, 2}
    pₛ::Real
end

function initialize_walker(N; y=N)
    x₀ = rand(1:N)
    y₀ = y
    return RandomWalker([Int8(y₀), Int8(x₀)])

end

function initialize_cluster(N, pₛ)
    all_idxs = 1:N 

    # set all indices to 1 initially
    #idx_dict = Dict{Tuple{Int8, Int8}, Int8}([(Int8(y), Int8(x)) => 1 for x in all_idxs, y in all_idxs])
    #seed_idx = (Int8(1), Int8(N/2))
    #idx_dict[seed_idx] = 0 # growth are characterized by 0

    domain = ones(Int8, N, N)
    domain[1, Int(N/2)] = 0 # position of initial seed
    return Cluster(N, domain, pₛ)
end

function walk!(walker::RandomWalker)
    choice = rand(keys(moves))     # choose random step
    walker.pos[:] += moves[choice] # update position
    return choice
end

function get_neighbours(walker, N)
    neighbours = ((walker.pos[1]+1, walker.pos[2]),
                  (walker.pos[1]-1, walker.pos[2]),
                  (walker.pos[1], walker.pos[2]+1),
                  (walker.pos[1], walker.pos[2]-1)) # neigbours of new cell

    # only choose neighbours who are not outside the domain
    return [Int8.(idx) for idx in neighbours if !(any(any.((idx .== 0, idx .== N+1))))]
end


function update!(walker, cluster)

    step = walk!(walker)
    where_cluster = findall(cluster.domain .== 0)
    if any(walker.pos[1] .== (0, cluster.N+1))
        walker.pos[:] = [N, rand(1:N)] # if walker exits domain, reset its position
    elseif walker.pos[2] == 0
        walker.pos[2] = N
    elseif walker.pos[2] == cluster.N+1
        walker.pos[2] = 1
    elseif CartesianIndex(w.pos...) in where_cluster
        walker.pos[:] -= moves[step]
    end


    neighbours = get_neighbours(walker, cluster.N)
    neighbour_cluster = any([iszero(cluster.domain[n...]) for n in neighbours])
    if neighbour_cluster
        if cluster.pₛ == 1 || cluster.pₛ >= rand()
            cluster.domain[CartesianIndex(walker.pos...)] = 0
            walker.pos[:] = [N, rand(1:N)] # remove walker from collection
        else    

            
        
        end
    end
end

function simulate(N, pₛ=1.0; num_steps=nothing, size_stop=nothing)

    ax = range(0, stop=1, length=N)

    cluster = initialize_cluster(N, pₛ)
    walker = initialize_walker(N)

    if !(isnothing(num_steps))
        for i = 1:num_steps
            update!(walker, cluster)
        end
    elseif !(isnothing(size_stop))
        clustersize = 1 - sum(cluster.domain)
        while true
            update!(walker, cluster)
            new_size = 1 - sum(cluster.domain)
            if new_size > clustersize
                clustersize = new_size
            end

            if clustersize >= size_stop
                break
            end
        end
    end

    return cluster
end


function animate_montecarlo_dla(num_steps; N=50, pₛ=1, step=1)
    ax = range(0, stop=1, length=N)

    cluster = initialize_cluster(N, pₛ)
    walker = initialize_walker(N)

    anim = @animate for i = 1:num_steps
        update!(walker, cluster)

        pos = (ax[walker.pos[2]], ax[walker.pos[1]])
        heatmap(ax, ax, cluster.domain)
        scatter!(pos, ylims=(0, 1), xlims=(0, 1), label=false)
    end every step

    gif(anim, "test_walker.gif")
end

# function montecarlo_dla(num_steps, nwalkers; N=50, pₛ=1, step=1)
#     walkers = Dict([i => initialize_walker(N) for i = 1:nwalkers])
#     ax = range(0, stop=1, length=N)

#     cluster, domain = initialize_domain(N)

#     anim = @animate for i = 1:num_steps
#         pos = Tuple{Float64,Float64}[]
#         for (id, walker) in walkers
#             walk!(walker)

#             if any(walker.pos[1] .== (0, N+1))
#                 #walker.pos[:] = [N, rand(1:N)]
#                 pop!(walkers, id)
#                 continue
#             elseif walker.pos[2] == 0
#                 walker.pos[2] = N
#             elseif walker.pos[2] == N+1
#                 walker.pos[2] = 1
#             end
#             #pos = walker.pos
#             try
#                 push!(pos, (ax[walker.pos[2]], ax[walker.pos[1]]))
#             catch e
#                 println(walker.pos)
#             end

#             neighbours = get_neighbours(walker, N)
#             neighbour_cluster = any([cluster.indices[n] == 0 for n in neighbours])
#             if neighbour_cluster
#                 if pₛ == 1 || pₛ >= rand()
#                     cluster.indices[Tuple(walker.pos)] = 0 # add to growth
#                     domain[CartesianIndex(walker.pos...)] = 0
#                     pop!(walkers, id) # remove walker from collection
#                 end
#             end
#             #println([ax[walker.pos[2]]], [ax[walker.pos[1]]])
#             #scatter!(p, [ax[walker.pos[2]]], [ax[walker.pos[1]]], label=false, color=:red)
#         end
#         heatmap(ax, ax, domain)
#         scatter!(pos, ylims=(0, 1), xlims=(0, 1), label=false)
#     end every step
#     gif(anim, "test_walker.gif")
# end