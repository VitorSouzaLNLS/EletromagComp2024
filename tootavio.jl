using Meshes
using GLMakie

esf = Sphere((0, 0, 0), 1)


fig = Figure()
ax = Axis3(fig[1,1], aspect=(1,1,1))
Meshes.viz!(ax, esf)
display(fig)
