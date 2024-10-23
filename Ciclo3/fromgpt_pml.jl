using LinearAlgebra, Plots

# Parameters
Nx, Ny = 70, 200  # Grid size
dx, dy = 0.1, 0.1  # Grid spacing
c = 1.0            # Wave propagation speed (normalized)
dt = 0.05          # Time step size
tmax = 600          # Total time steps

# Source: position and frequency of the wave source
src_x, src_y = Nx ÷ 3, Ny ÷ 3  # Source position at center
freq = 1.0

# Fields: initialize electric and magnetic fields
E = zeros(Float64, Nx, Ny)  # Electric field (E-field)
Hx = zeros(Float64, Nx, Ny) # Magnetic field component (x-direction)
Hy = zeros(Float64, Nx, Ny) # Magnetic field component (y-direction)

if true
    pml = ones(Nx, Ny)
    pml_xwidth = 7
    pml_ywidth = 10
    K = 2
    pml_stren = 1.5
    for i in 1:Nx
        dist = 0.0
        if i<pml_xwidth
            dist = pml_xwidth-i
        elseif i>Nx-pml_xwidth
            dist = i - (Nx - pml_xwidth) -1
        end
        # pml[i,:] *= 1 - (dist/pml_xwidth)^K
        pml[i,:] *= exp(-pml_stren*(dist/pml_xwidth)^2)
    end

    for j in 1:Ny
        dist = 0.0
        if j < pml_ywidth
            dist = pml_ywidth - j
        elseif j>Ny-pml_ywidth
            dist = j - (Ny - pml_ywidth) -1
        end
        # pml[:,j] *= 1 - (dist/pml_ywidth)^K
        pml[:,j] *= exp(-pml_stren*(dist/pml_xwidth)^2)
    end

    heatmap(pml)
end

if true
    nref = ones(Float64, Nx, Ny)
    prism_pos = src_y + (Ny ÷ 3)
    prism_width = 20
    for j in prism_pos-prism_width:prism_pos+prism_width
        if j>=prism_pos
            xmax = Nx - Int(ceil((j + 1 - prism_pos)*Nx/prism_width))
        else
            xmax = 2 + Int(1 + ceil((j-prism_pos+prism_width)*(Nx/prism_width)))
        end
        for i in 1:xmax
            nref[i, j] = 0.5
        end
    end
    # nref[src_x, src_y] = 0
    heatmap(nref)
end

# Update the electric and magnetic fields
function update_fields!(E, Hx, Hy)
    # Update the magnetic field (Hx, Hy) based on the curl of the E-field

    for i in 1:Nx-1, j in 1:Ny-1
        Hx[i, j] += pml[i, j] * dt * (E[i, j+1] - E[i, j]) / dy
        Hy[i, j] += pml[i, j] * dt * (E[i+1, j] - E[i, j]) / dx
    end

    # Update the electric field (E) based on the curl of the magnetic field
    for i in 2:Nx, j in 2:Ny
        E[i, j] += dt * c*nref[i,j] * ((Hy[i, j] - Hy[i-1, j]) / dx + (Hx[i, j] - Hx[i, j-1]) / dy)
        # Apply PML damping
        E[i, j] *= pml[i, j]
    end
end

# Main loop: propagate the wave

# E[src_x, src_y] = 1.0
# freq = 3.0

osci = true

if !osci
    for i in 1:Nx
        for j in 1:Ny
            E[i, j] = 3 * exp(-((src_x - i)^2 + (src_y - j)^2)/10)
        end
    end
end

anim = @animate for t in 1:tmax
    # Add source: oscillating electric field at the center
    # if t < 4
    if osci
        E[src_x, src_y] = sin(2π * freq * t * dt)
    end

    # Update fields
    update_fields!(E, Hx, Hy)

    # Plot the electric field at each time step
    heatmap(E, clim=(-1, 1), color=:inferno, size=(Ny*5, Nx*5))
    plot!([prism_pos-prism_width, prism_pos, prism_pos+prism_width, prism_pos-prism_width], [0, Nx, 0, 0], color=:black, legend=false, fill=(0, 0.3, :gray))
end

# Display
gif(anim, "em_wave_simulation_pml.gif", fps=1/dt)
