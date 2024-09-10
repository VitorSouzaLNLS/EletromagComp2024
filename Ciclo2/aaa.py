import numpy as np
import matplotlib.pyplot as plt

plt.switch_backend('TKAgg')

# Example function that returns a vector (vx, vy, vz) for a given point (x, y, z)
def vector_field(x, y, z):
    # Simple example: radial field pointing away from origin
    vx = x
    vy = y
    vz = z
    return vx, vy, vz

# Generate a grid of points (x, y, z) where you want to plot the vectors
points = [(x, y, z) for x in np.linspace(-1, 1, 5) for y in np.linspace(-1, 1, 5) for z in np.linspace(-1, 1, 5)]

# Unpack the points into x, y, z lists
x, y, z = zip(*points)

# Calculate the vectors at each point
vx, vy, vz = zip(*[vector_field(xi, yi, zi) for xi, yi, zi in points])

# Create 3D figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the vectors using quiver
ax.quiver(x, y, z, vx, vy, vz, length=0.2, normalize=True)

# Label axes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
