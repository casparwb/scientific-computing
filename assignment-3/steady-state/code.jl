
using Plots, PlotThemes, LinearAlgebra, SparseArrays
using Arpack: eigs
theme(:juno)

function LaplacianMatrix(Nx, Ny, Lx, Ly; return_sparse=true)
    dx = Lx / Nx
    dy = Ly / Ny

    #Dx = Tridiagonal(-1ones(Nx-1), ones(Nx), zeros(Nx-1)) / dx
    #Dy = Tridiagonal(-1ones(Ny-1), ones(Ny), zeros(Ny-1)) / dy

    Dx = [ [1.0 zeros(1,Nx-1)]; diagm(1 => ones(Nx-1)) - I(Nx)] / dx
    Dy = [ [1.0 zeros(1,Ny-1)]; diagm(1 => ones(Ny-1)) - I(Ny)] / dy


    Ax = Dx' * Dx
    Ay = Dy' * Dy

    if return_sparse
        return -1*sparse(kron(I(Ny), Ax) + kron(Ay, I(Nx)))
    else
        return -1*(kron(I(Ny), Ax) + kron(Ay, I(Nx)))
    end
 end


 function circle(Nx, Ny, D; sparse_matrix=false)

    R = D/2

    x = collect(range(-R, stop=R, length=Nx))
    y = collect(range(-R, stop=R, length=Ny))'

    r = sqrt.(x .^ 2 .+ y .^ 2)
    
    """circle_idxs = findall(r .< R)
    rows = [idx[1] for idx in circle_idxs]
    cols = [idx[2] for idx in circle_idxs]"""

    i = [i for i in eachindex(r) if r[i] < R]

    M = LaplacianMatrix(Nx, Ny, D, D, return_sparse=sparse_matrix)

    M = M[i, i]
    if sparse_matrix
        e = eigs(M, which=:SM)
    else
        e = eigen(M)
    end

    return e, i
    #return r
end

function steady_state_concentration(Nx, Ny; radius=2, source=(0.6, 1.2))

    c = zeros(Nx, Ny)

    x_ax = range(-radius, stop=radius, length=Nx)
    y_ax = range(-radius, stop=radius, length=Ny)


    source_x_idx = argmin(abs.(x_ax .- source[1]))
    source_y_idx = argmin(abs.(y_ax .- source[2]))

    c[source_y_idx, source_x_idx] = 1.0
    r = sqrt.(x_ax .^ 2 .+ y_ax' .^ 2)

    M = LaplacianMatrix(Nx, Ny, 2*radius, 2*radius)

    i = [i for i in eachindex(r) if r[i] < radius]

    M = M[i, i]

    b = reshape(c, :, 1) # b is c collapsed column-wise
    # b = zeros(size(M, 2))

    x = M \ -b[i]

    c = zeros(Nx, Ny)
    c[i] = x

    return x_ax, y_ax, c
end

function plot_steady_state(x, y, c; radius=2)


    θ = range(0, stop=2π, length=100)

    heatmap(x, y, c, aspect_ratio=1, xlabel="x", ylabel="y", colorbar_title="Concentration")
    
    plot!(radius*cos.(θ), radius*sin.(θ), label="Domain")

end