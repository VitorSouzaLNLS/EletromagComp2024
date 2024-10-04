using GLMakie, Meshes
using FileIO, MeshIO, MeshBridge

# Create the sphere
esf = Meshes.Sphere((0.0, 0.0, 0.0), 1.0)
catload = load("12222_Cat_v1_l3.obj")
catmesh = convert(Meshes.Mesh, catload)


# Create the figure and axis
if true
    fig = Figure()
    ax = Axis3(fig[1, 1])
    # sc = Scene()
    Meshes.viz!(ax, catmesh)
    Meshes.viz!(ax, esf)

    # Rotate and save each frame
    n_frames = 30  # Number of frames
    fps = 10        # Frames per second
    # cam = GLMakie.Camera3D(sc)

    # Loop through each frame
    # for i in 1:n_frames
    #     # Rotate the camera around the z-axis
    #     # ax.scene.camera.rotation[] = Makie.Vec3f0(0, 0, 2π * i / n_frames)


    #     # Save each frame as an image
    d = 100
    xlims!(ax, -d/2, d/2)
    ylims!(ax, -d/2, d/2)
    zlims!(ax, -d/5, 6*d/5)

    start_angle = π / 4
    n_frames = 120
    ax.viewmode = :fit # Prevent axis from resizing during animation
    record(fig, "testgif/test.mp4", 1:n_frames) do frame
        println(frame)
        ax.azimuth[] = start_angle + 2pi * frame / n_frames
    end
end

# using GLMakie

# xs = -10:0.5:10
# ys = -10:0.5:10
# zs = [cos(x) * sin(y) for x in xs, y in ys]

# fig = GLMakie.Figure()
# ax = GLMakie.Axis3(fig[1, 1])
# Makie.surface!(xs, ys, zs)

# start_angle = π / 4
# n_frames = 120
# ax.viewmode = :fit # Prevent axis from resizing during animation
# record(fig, "test.mp4", 1:n_frames) do frame
#     ax.azimuth[] = start_angle + 2pi * frame / n_frames
# end
