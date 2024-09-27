using MeshIO
using FileIO
using GLMakie

omesh = load("Ciclo2/airplane_model/11805_airplane_v2_L2.obj")


fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1.7,2.2,1))
mesh!(ax, omesh)
display(fig)
