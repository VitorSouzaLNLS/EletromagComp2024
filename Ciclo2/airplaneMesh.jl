# using GeometryBasics
using MeshIO
using FileIO
using GLMakie: Figure, mesh!, Axis3, GeometryBasics, quiver!, xlims!, zlims!, ylims!, arrows!, Vec3, lines!, scatter!
using Meshes
using MeshBridge

function is_outside_sphere(x_::Float64, y_::Float64, z_::Float64, radius::Float64)
    r_ = sqrt(x_^2 + y_^2 + z_^2)  # Distância radial no plano
    return r_ > radius  # Retorna true se o ponto estiver fora do cone
end

struct element
    center
    area
end

mutable struct surface
    elements::Vector{element}
    Q::Float64
    surface(elems::Vector{element}, Q=0.0) = new(elems, Q)
end

function surface(m::Mesh, Q=0.0)
    elems = element[]
    for f in faces(m, 2)
        push!(elems, element(coordinates(centroid(f)), area(f)))
    end
    return surface(elems, Q)
end

get_area(s::surface) = length(s.elements) == 0 ? 0.0 : sum([e.area for e in s.elements])

struct charge
    r::Tuple{Float64, Float64, Float64}
    Q::Float64
end

norm(vec) = sqrt(sum(vec.^2))

function calc_E_field(s::surface, r::Tuple{Float64, Float64, Float64}, charges::Vector{charge}=charge[])
    E = [0.0, 0.0, 0.0]
    tol = 0.0
    A = get_area(s)
    # println(A)
    for e in s.elements
        rvec = [r[1] - e.center[1], r[2] - e.center[2], r[3] - e.center[3]]
        q = (e.area / A) * s.Q
        ramp = norm(rvec)^3
        if ramp > tol
            E = E .+ (q .* rvec ./ ramp)
        end
    end
    for e in charges
        rvec = [r[1] - e.r[1], r[2] - e.r[2], r[3] - e.r[3]]
        ramp = norm(rvec)^3
        if ramp > tol
            E = E .+ (e.Q .* rvec ./ ramp)
        end
    end
    return E
end

# omesh = load("Ciclo2/airplane_model/11805_airplane_v2_L2.obj")
# m = convert(Mesh, omesh)

R = 1.0
c = Meshes.Point((0.0, 0.0, 0.0))
s = Meshes.Sphere(c, R)
m = discretize(s, RegularDiscretization((20, 20)))

# surf = surface(m, 1.0)
surf = surface(element[], 0.0)
cargas = [
    charge((-1.0, 0.0, 0.0), +1.0),
    # charge((1.0, 0.0, 0.0), -1.0),
    charge((0.0, 1.0, 0.0), -1.0)
]


# x_range = range(-1.5, 1.5, length=9)
# y_range = range(-1.5, 1.5, length=9)
# z_range = range(-1.0, 1.0, length=3)
# gr = [(x, y, z) for x in x_range for y in y_range for z in z_range]
# vx, vy, vz, x, y, z, stren = [], [], [], [], [], [], []

# for (xi, yi, zi) in gr
#     if is_outside_sphere(xi, yi, zi, R)
#         E = calc_E_field(surf, (xi, yi, zi), cargas)
#         push!(vx, E[1])
#         push!(vy, E[2])
#         push!(vz, E[3])
#         push!(x, xi)
#         push!(y, yi)
#         push!(z, zi)
#         push!(stren, norm(E))
#     end
# end

# calculate 1 line
function calc_line(s::surface, charges::Vector{charge}, r0::Tuple)
    points = []
    push!(points, Vector([r0...]))
    ds = 0.1
    thr = 0
    while thr < 2
        p = points[end]
        E = calc_E_field(s, (p[1], p[2], p[3]), charges)
        # println(norm(E ./ norm(E)))
        push!(points, p .+ ((E ./ norm(E)) .* ds))
        thr += 1#norm(points[end] .- points[end-1])
    end
    return points
end

starts = [
    # (-0.8, 0.0, 0.0),
    # (-1.0, 0.0, +0.2),
    # (-1.0, 0.0, -0.2),
    # (-1.0, -0.2, 0.0),
    (-0.2, +0.2, 0.6)
]

if true
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=(1,1,1))
    # mesh!(ax, GeometryBasics.Sphere(GeometryBasics.Point3((1.0, 2.0, 3.0))))
    # quiver!(x, y, z, vx, vy, vz)
    # viz!(m)
    for q in cargas
        viz!(Meshes.Sphere(Meshes.Point3(q.r...), 0.2))
    end

    # arrows!(x, y, z, vx, vy, vz, lengthscale=0.3, arrowsize=Vec3((0.1, 0.1, 0.1)))

    for start in starts
        x, y, z = [], [], []
        pf = calc_line(surf, cargas, start)
        # x, y, z = [], [], []
        for p in pf
            push!(x, p[1])
            push!(y, p[2])
            push!(z, p[3])
        end
        println(length(x))
        scatter!(ax, x, y, z)
    end

    xlims!(ax, -2, 2)
    ylims!(ax, -2, 2)
    zlims!(ax, -2, 2)

    display(fig)
end
