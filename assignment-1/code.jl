### Set 1: Vibrating String ###

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


function run(string, T)

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

function animate(string, T)

    t = 0.0
    N = Int(T/string.Δt)
    amplt = maximum(string.Ψ[:,1])
    anim = @animate for i ∈ 1:N
        step!(string)
        t += string.Δt
        plot(string.Ψ[:,3], ylims=(-amplt, amplt))
    end

    gif(anim, "test_string.gif", fps = 10)
end