using PyCall

mplt = pyimport("matplotlib")
plt = pyimport("matplotlib.pyplot")
np = pyimport("numpy")


mplt.use("tkagg")

function plot_cone(R, H, ax, subdiv=10)
    # Make data
    N = Int(ceil(subdiv))
    h = Float64(H)
    theta = np.linspace(0, 2π, N)  # Angular divisions
    z = np.linspace(0, h, N)       # Height divisions

    # Create arrays for x, y, z coordinates
    x = np.zeros((length(z), length(theta)))
    y = np.zeros((length(z), length(theta)))
    Z = np.zeros((length(z), length(theta)))

    for i in eachindex(z)
        r = R * (1 - z[i] / h)  # Radius decreases linearly with height
        for j in eachindex(theta)
            x[i, j] = r * cos(theta[j])
            y[i, j] = r * sin(theta[j])
            Z[i, j] = z[i]
        end
    end

    # Plot the surface
    ax.plot_surface(x, y, Z)
end


# Call the function with desired base radius R and height H
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
plot_cone(5, 10, ax)
# ax.plot([0.0,], [0.0,], [0.0,], "o")
ax.set_aspect("equal")
plt.show()
