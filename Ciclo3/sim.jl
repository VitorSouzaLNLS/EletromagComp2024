using Plots

# Parameters
Nx = 100
Ny = 50
c = 1.0
Δx = 1.0
Δy = Δx
Δt = 0.5
T = Int(ceil((Nx) / Δt))

nindex = ones(Nx, Ny)
nindex[60:end, :] .= 0.5

# Initialize fields
u = zeros(Nx, Ny)
u_prev = zeros(Nx, Ny)
u_next = zeros(Nx, Ny)

# Initial conditions (Gaussian pulse in the center)
x0, y0 = Nx ÷ 2, Ny ÷ 2
for i in 1:Nx
    for j in 1:Ny
        u[i, j] = 3*exp(-((i - x0)^2 + (j - y0)^2) / 10.0)
    end
end

u_prev .= u

plt = plot(size=(600*(Nx/100), 600*(Ny/100)))
pml = ones(Nx, Ny)
pml_xwidth = 10
pml_ywidth = 10
K = 2.0
for i in 1:Nx
    dist = 0.0
    if i<=pml_xwidth
        dist = pml_xwidth-i + 1
    elseif i>=Nx-pml_xwidth
        dist = pml_xwidth - (Nx - i)
    end
    pml[i,:] *= 1 - (dist/pml_xwidth)^K
    # pml[i,:] *= exp(-(dist^2)/pml_width / 2)
end

for j in 1:Ny
    dist = 0.0
    if j<=pml_ywidth
        dist = pml_ywidth-j + 1
    elseif j>=Ny-pml_ywidth
        dist = pml_ywidth - (Ny - j)
    end
    pml[:,j] *= 1 - (dist/pml_ywidth)^K
end

# Time stepping with PML
for t in 1:T
    for i in 2:Nx-1
        for j in 2:Ny-1
            u_next[i, j] = (2 - Δt * (1 - pml[i, j]))*u[i, j] - (1 - Δt * (1 - pml[i, j]))*u_prev[i, j] +
            ((c*nindex[i,j])^2 * (Δt*pml[i,j])^2 / Δx^2) * (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) +
            ((c*nindex[i,j])^2 * (Δt*pml[i,j])^2 / Δy^2) * (u[i, j+1] - 2 * u[i, j] + u[i, j-1])
        end
    end

    # Update the fields for the next iteration
    u_prev .= u
    u .= u_next

    # Plot every few time steps
    ucp = transpose((zeros(Nx, Ny) + u)[pml_xwidth:end-pml_xwidth, pml_ywidth:end-pml_ywidth])
    if t % 2 == 0
        heatmap!(ucp, title="Time t = $(t*Δt)", color=:inferno, clims=(-1, 1))
        display(plt)
        # sleep(0.05)
    end
end
