using LinearAlgebra, Plots

# Parameters
Nx, Ny = 100, 100  # Grid size
dx, dy = 0.1, 0.1  # Grid spacing
c = 1.0            # Wave propagation speed (normalized)
dt = 0.05          # Time step size
tmax = 150          # Total time steps

# Fields: initialize electric and magnetic fields
E = zeros(Float64, Nx, Ny)  # Electric field (E-field)
Hx = zeros(Float64, Nx, Ny) # Magnetic field component (x-direction)
Hy = zeros(Float64, Nx, Ny) # Magnetic field component (y-direction)

# Source: position and frequency of the wave source
src_x, src_y = Nx ÷ 2, Ny ÷ 2  # Source position at center
freq = 1.0                     # Source frequency

# Update the electric and magnetic fields
function update_fields!(E, Hx, Hy)
    # Update the magnetic field (Hx, Hy) based on the curl of the E-field
    for i in 1:Nx-1, j in 1:Ny-1
        Hx[i, j] += dt * (E[i, j+1] - E[i, j]) / dy
        Hy[i, j] += dt * (E[i+1, j] - E[i, j]) / dx
    end

    # Update the electric field (E) based on the curl of the magnetic field
    for i in 2:Nx, j in 2:Ny
        E[i, j] += dt * c * ((Hy[i, j] - Hy[i-1, j]) / dx + (Hx[i, j] - Hx[i, j-1]) / dy)
    end
end

# Main loop: propagate the wave
anim = @animate for t in 1:tmax
    # Add source: oscillating electric field at the center
    E[src_x, src_y] = sin(2π * freq * t * dt)

    # Update fields
    update_fields!(E, Hx, Hy)

    # Plot the electric field at each time step
    heatmap(E, clim=(-1, 1), color=:blues, title="Electric Field at time step $t")
end

# Display the animation
gif(anim, "em_wave_simulation.gif", fps=20)
