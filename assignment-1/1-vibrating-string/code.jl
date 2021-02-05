### Set 1: Vibrating String ###
using Plots, DifferentialEquations

struct StringConstruct
    c::Float64
    Δx::Float64
    Δt::Float64
    Γ::Float64
    t::Array{Float64}
    Ψ::Array{Float64, 2}
    boundaries::Tuple{Float64, Float64}
end

function set_boundaries!(string::StringConstruct, boundaries)
    string.Ψ[[1, end], :] .= boundaries
end

function initialize(init_func, N, Δt, c=1.0, L=1.0, boundaries=(0.0, 0.0))

    Δx = L/N
    x = range(0, stop=L, length=N)

    Ψ = zeros(Float64, N, 3)
    Ψ[:,1] = init_func(x)
    Ψ[:,2] = Ψ[:,1]

    Γ = (Δt*c/Δx)^2
    string = StringConstruct(c, Δx, Δt, Γ, [0.0], Ψ, boundaries)

    set_boundaries!(string, boundaries)
    return string
end


function step!(string)

    term1 = circshift(string.Ψ[:,2], -1) - 2*string.Ψ[:,2] + circshift(string.Ψ[:,2], 1)
    term2 = 2*string.Ψ[:,2] - string.Ψ[:,1]
    string.Ψ[:,3] = string.Γ*term1 + term2
    set_boundaries!(string, string.boundaries)

    string.Ψ[:,1] = string.Ψ[:,2]
    string.Ψ[:,2] = string.Ψ[:,3]
end


function simulate(string, T)

    Ψvals = []
    push!(Ψvals, string.Ψ[:,2])
    t = 0.0
    while t < T
        step!(string)
        push!(Ψvals, string.Ψ[:,3])
        t += string.Δt
    end

    return Ψvals
end


function animate(string::StringConstruct, T, outpath, fps=10, saveevery=10)

    N = Int(T/string.Δt)
    t = 0.0
    anim = @animate for i ∈ 1:N
        step!(string)
        t += string.Δt
        plot(string.Ψ[:,3], ylims=(-1, 1), title="t=$(round(t, digits=2))", legend=false)
        xlabel!("x")
        ylabel!("Ψ")
    end

    gif(anim, outpath*".gif", fps = fps)
end

function animate(Ψ, T, outpath; fps=10)
    t = range(0, stop=T, length=size(Ψ, 1))
    x = range(0, stop=1, length=size(Ψ[1], 1))
    anim = @animate for i ∈ 1:length(t)
        plot(x, Ψ[i], ylims=(-1, 1), title="t=$(round(t[i], digits=2))", legend=false)
        xlabel!("x")
        ylabel!("Ψ")
    end

    gif(anim, outpath*".gif", fps = fps)
end

function plot_wave(Ψ, T, outpath, title; nplots=10)

    t = range(0, stop=T, length=size(Ψ, 1))
    x = range(0, stop=1, length=size(Ψ[1], 1))

    step = round(Int, size(Ψ, 1)/nplots)
    #tlabels = ["t = $(round(t[i], digits=2))" for i = 1:step:length(t)]
    p = plot()
    for i ∈ 1:step:length(t)
        label = "t = $(round(t[i], digits=2))"
        plot!(p, x, Ψ[i], ylims=(-1, 1), label=label)
    end
    xlabel!("x")
    ylabel!("Ψ")
    title!(title)
    savefig(outpath*".svg")
end

function visualize(N, Δt, T, fps=10, nplots=10)

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

    results = []
    for (init, outname) in zip(init_funcs, outnames)
        string = initialize(init, N, Δt)
        Ψ = simulate(string, T)
        #animate(Ψ, T, outname)
        plot_wave(Ψ, 0.2, outname, "outname", nplots=5)
    end


end
