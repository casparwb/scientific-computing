using Plots, PlotThemes, LinearAlgebra, SparseArrays, BenchmarkTools
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

function InitializeWave(Nx, Ny)

    domain = zeros(Float64, Ny, Nx)

    domain[round(Int, Ny/2), round(Int, Nx/2)] = 1

    return reshape(domain, :, 1)

end

function EigProb(Nx, Ny, Lx, Ly; sparse_matrix=false)

    M = LaplacianMatrix(Nx, Ny, Lx, Ly, return_sparse=sparse_matrix)

    if sparse_matrix
        eig = eigs(M, which=:SM)
    else
        eig = eigen(M)
    end

    return eig
end

function plot_eigvecs(eigvec, Nx, Ny, Lx, Ly)

    x = range(0, stop=Lx, length=Nx)
    y = range(0, stop=Ly, length=Ny)

    vec_reshaped = reshape(eigvec, Ny, Nx)

    heatmap(x, y, vec_reshaped)

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

function plot_circle_eigvec(Nx, Ny, L, eig, circle_idxs, mode)



    x_ax = range(-L/2, stop=L/2, length=Nx)
    y_ax = range(-L/2, stop=L/2, length=Ny)

    eigvec = zeros(Nx, Ny)

    eigvec[circle_idxs] = eig[2][:,mode]
    #eigvec[circle_idxs] = eig.vectors[:,mode]

    heatmap(x_ax, y_ax, eigvec, aspect_ratio=1)

end

function K_to_λ(K)
    return sqrt(-K)
end

function plot_eigvec_shapes(Nx = 100, Ny = 100)

    L = 1.0

    ### a square ###
    Lx = Ly = L
    ax = range(0, stop=L, length=Nx)
    eig = EigProb(Nx, Ny, Lx, Ly, sparse_matrix=true)

    eigvec1 = eig[2][:,1]
    λ₁ = K_to_λ(eig[1][2])
    eigvec2 = eig[2][:,2]
    λ₂ = K_to_λ(eig[1][3])

    p1 = heatmap(ax, ax, reshape(eigvec1, Nx, Ny), title="Square | λ₂=$(round(λ₁, digits=3))")
    p2 = heatmap(ax, ax, reshape(eigvec2, Nx, Ny), title="Square | λ₃=$(round(λ₂, digits=3))")

    ### rectangle ###
    Lx = 2*L
    Ly = L

    x_ax = range(0, stop=Lx, length=Nx)
    y_ax = range(0, stop=Ly, length=Ny)
    eig = EigProb(Nx, Ny, Lx, Ly, sparse_matrix=true)

    eigvec1 = eig[2][:,2]
    λ₁ = K_to_λ(eig[1][2])
    eigvec2 = eig[2][:,3]
    λ₂ = K_to_λ(eig[1][3])

    p3 = heatmap(x_ax, y_ax, reshape(eigvec1, Nx, Ny), title="Rectangle | λ₂=$(round(λ₁, digits=3))")
    p4 = heatmap(x_ax, y_ax, reshape(eigvec2, Nx, Ny), title="Rectangle | λ₃=$(round(λ₂, digits=3))")

    ### circle ###
    Lx = Ly = L

    x_ax = range(-Lx/2, stop=Lx/2, length=Nx)
    y_ax = range(-Ly/2, stop=Ly/2, length=Ny)

    eig, circle_idxs = circle(Nx, Ny, L, sparse_matrix=true)

    eigvec1 = zeros(Nx, Ny)
    eigvec2 = zeros(Nx, Ny)

    eigvec1[circle_idxs] = eig[2][:,1]
    λ₁ = K_to_λ(eig[1][1])
    eigvec2[circle_idxs] = eig[2][:,3]
    λ₂ = K_to_λ(eig[1][3])

    p5 = heatmap(x_ax, y_ax, eigvec1, title="Circle | λ₁=$(round(λ₁, digits=3))", aspect_ratio=1)
    p6 = heatmap(x_ax, y_ax, eigvec2, title="Circle | λ₃=$(round(λ₂, digits=3))", aspect_ratio=1)

    plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(1000, 900))
end

function benchmark_eig(Nx = 100, Ny = 100; verbose=true)

    Lx = Ly = 1

    function output(bm)

        println("Minimum time: $(minimum(bm))")
        println("Median time: $(median(bm))")
        println("Mean time: $(mean(bm))")
        println("Maximum time: $(maximum(bm))")

    end

    ### benchmark non-sparse ###

    benchmark_non_sparse = @benchmark EigProb($Nx, $Ny, $Lx, $Ly, sparse_matrix=false)

    ### benchmark sparse ###

    benchmark_sparse = @benchmark EigProb($Nx, $Ny, $Lx, $Ly, sparse_matrix=true)

    if verbose
        println("Non-sparse result:")
        output(benchmark_non_sparse)
        
        println("\nSparse result:")
        output(benchmark_sparse)
    end

    return benchmark_non_sparse, benchmark_sparse

end

function eigenfrequency_spectrum(Nx=100, Ny=100)

    L_values = 1:10
    λs = zeros(3, length(L_values))

    for (i, L) in enumerate(L_values)
        ### square ###
        Lx = Ly = L
        eig = EigProb(Nx, Ny, Lx, Ly, sparse_matrix=true)
    
        λ = K_to_λ(eig[1][1])

        λs[1,i] = λ

        ### rectangle ###
        Lx = 2*L
        Ly = L
        eig = EigProb(Nx, Ny, Lx, Ly, sparse_matrix=true)
        λ = K_to_λ(eig[1][1])
        λs[2,i] = λ

        ### circle ###
        Lx = Ly = L
        eig = circle(Nx, Ny, L, sparse_matrix=true)[1]
        λ = K_to_λ(eig[1][1])
        λs[3,i] = λ
    end

    return L_values, λs
end

function time_dependent_solution(Nx = 100, Ny = 100; eigidx=1)

    L = 1
    Lx = 2*L
    Ly = 1

    #eig = EigProb(Nx, Ny, Lx, Ly, sparse_matrix=true)
    eig, idxs = circle(Nx, Ny, L, sparse_matrix=true)

    K = eig[1]
    v = eig[2]

    T(t, λ) = cos(λ*t) + sin(λ*t)

    t_values = range(0, stop=10, length=1000)

    λ = K_to_λ(K[eigidx])
    #eigvec = reshape(v[:,eigidx], Nx, Ny)

    eigvec = zeros(Nx, Ny)

    eigvec[idxs] = v[:,eigidx]
    
    anim = @animate for t in t_values

        u = eigvec*T(t, λ)

        surface(u, colorbar=false, zlims=(-0.05, 0.05))
    end

    gif(anim, "circle_mode_$(eigidx).gif")
end