function finite_difference(Ψ, Γ)
    diff = circshift(Ψ, -1) - 2*Ψ + circshift(Ψ, 1)
    diff[[1, end]] .= 0.0
    return Γ*diff
end

function wave_integ(du, u, p, t)
    Γ = p
    ∂²Ψ∂x² = finite_difference(u[1], Γ)
    du[1] = u[2]
    du[2] = ∂²Ψ∂x²
end


function solve_wave(f, L, c, N, T)
    Δx = L/N
    x = range(0, stop=L, length=N)

    u0 = f(x)
    u0[[1, end]] .= 0.0
    u1 = zeros(Float64, N)
    Γ = (c/Δx)^2


end