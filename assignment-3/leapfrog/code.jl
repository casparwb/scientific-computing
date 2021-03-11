using Plots, PlotThemes
theme(:juno)


function F(x, t, k)
    return -k*x
end

function LeapfrogIntegration(N, k, Δt; x₀=1, v₀=0, f=F)


    x₀ = 1
    x = zeros(N)
    v = zeros(N)
    t = zeros(N)
    
    v[1] = v₀ .+ 0.5*F(x₀, k)
    x[1] = x₀

    for i = 1:N-1
        t[i+1] = (i+1)*Δt
        a₀ = f(x[i], t[i], k)

        x[i+1] = x[i] + v[i]*Δt + 0.5*a₀*Δt^2

        a₁ = f(x[i+1], t[i+1], k)

        v[i+1] = v[i] + 0.5*(a₀ + a₁)*Δt # mean of old and new acceleration times time step


    end

    return x, v

end

function F_RK4(u, k)
    return [u[2], F(u[1], k)]
end

function RK45(N, k, Δt; x₀=1, v₀=0)

    u = zeros(2, N)
    
    v₁ = v₀ + 0.5*F(x₀, k)

    u[:, 1] = [x₀, v₁]


    for i = 1:N-1
        K1 = Δt*F_RK4(u[:, i],           k)
        K2 = Δt*F_RK4(u[:, i] .+ 0.5*K1, k)
        K3 = Δt*F_RK4(u[:, i] .+ 0.5*K2, k)
        K4 = Δt*F_RK4(u[:, i] .+ K3,     k)

        u[:, i+1] = u[:, i] .+ 1/6*(K1 .+ K2 .+ K3 .+ K4)
    end

    return u[1,:], u[2,:]
end

function vary_k(N, Δt; x₀=1, v₀=0)


    k_values = range(0.1, stop=1, length=5)

    x_results = zeros(length(k_values), N)
    v_results = similar(x_results)

    for (i, k) in enumerate(k_values)
        t, x, v = LeapfrogIntegration(N, k, Δt, x₀ = x₀, v₀ = v₀)

        x_results[i,:] = x
        v_results[i,:] = v

    end

    t = range(0, stop=Δt*N, length=N)

    labels = reshape(["k = $(round(k, digits=2))" for k in k_values], 1, :)
    p1 = plot(t, x_results', ylabel="x", label=labels)

    p2 = plot(t, v_results', xlabel="t", ylabel="v", label=false)

    plot(p1, p2, layout=(2, 1), size=(800, 700))

end

function compare_solvers(N, k, Δt; x₀=1, v₀=0)

    results_leapfrog = LeapfrogIntegration(N, k, Δt, x₀=x₀, v₀=v₀)[1]
    results_rk4 = RK45(N, k, Δt, x₀=x₀, v₀=v₀)[1]


    t = range(0, stop=N*Δt, length=N)

    p1 = plot(t, results_leapfrog, label="Leapfrog")
    plot!(t, results_rk4, label="RK4")
    xlabel!("t")

end

    

function F_vel(x, t, k, ω)
    return sin(ω*t) - k*x
end

function sinusoidal_force(N, k, Δt; x₀=1, v₀=0)

    true_ω = sqrt(k)
    ωs = collect(range(true_ω/2, stop=true_ω*5, length=5))
    push!(ωs, true_ω)


    results = zeros(2, length(ωs), N)

    p = plot()
    t = range(0, stop=N*Δt, length=N)
    for (i, ω) in enumerate(ωs)
        ff(x, t, k) = F_vel(x, t, k, ω)
        results = LeapfrogIntegration(N, k, Δt, f=ff)

        if i == length(ωs)
            label = "true ω = $(round(true_ω, digits=2))"
        else
            label="ω = $(round(ω, digits=2))"
        end

        plot!(p, results[1], results[2], label=label, xlabel="x", ylabel="v")
    end

    # ff(x, t, k) = F_vel(x, t, k, true_ω)
    # results = LeapfrogIntegration(N, k, Δt, f=ff)
    # plot!(p, results[1], results[2], label="true ω = $(round(true_ω, digits=2))", xlabel="x", ylabel="v")


    return p
end

