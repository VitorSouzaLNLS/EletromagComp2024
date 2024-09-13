using PyCall

mplt = pyimport("matplotlib")
plt = pyimport("matplotlib.pyplot")
cm = pyimport("matplotlib.cm")
np = pyimport("numpy")

mplt.use("tkagg")

struct element
    center
    area
end

mutable struct surface
    elements::Vector{element}
    Q::Float64
    surface(elems, Q=0.0) = new(elems, Q)
end

area(s::surface) = length(s.elements) > 0 ? sum([e.area for e in s.elements]) : 0.0

function frustum_area(R, r, dz)
    L = sqrt(dz^2 + (R - r)^2)
    A = pi * L * (R + r)
    return A
end

function cone_area(r, h)
    L = sqrt(h^2 + r^2)
    return pi * r * L
end

circle_area(radius) = pi * radius^2

sector_area(R, r) = pi * (R^2 - r^2)

function gen_cone(radius, height, num_segs)
    elements = element[]
    N = num_segs
    theta = [2pi*j/N for j in 0:N-1]

    # lateral
    dz = height / N * 2
    z = range(start=0, stop=height, step=dz)
    r = radius .* (1.0 .- (z ./ height))
    for i in 1:length(z)-1
        Ri, ri = r[i], r[i+1]
        a = frustum_area(Ri, ri, dz) ./ length(theta)
        zi = z[i] + dz/2
        x = (Ri + ri)/2 .* cos.(theta)
        y = (Ri + ri)/2 .* sin.(theta)
        for (xi, yi) in zip(x,y)
            push!(elements, element((xi, yi, zi), a))
        end
    end

    # base
    rings = Int(ceil(N/2))
    dr = radius/rings
    r = [radius-(i*dr) for i in 0:rings]
    for i in 1:length(r)-1
        Ri, ri = r[i], r[i+1]
        if ri == 0
            a = circle_area(Ri - ri)
            push!(elements, element((0,0,0), a))
        else
            a = sector_area(Ri, ri)/N
            x = (Ri + ri)/2 .* cos.(theta)
            y = (Ri + ri)/2 .* sin.(theta)
            for (xi, yi) in zip(x,y)
                push!(elements, element((xi, yi, 0), a))
            end
        end
    end

    return surface(elements)
end

function is_outside_cone(x_::Float64, y_::Float64, z_::Float64, radius::Float64, height::Float64)
    if (z_ <= 0.0) || (z_ >= height)
        return true  # Ponto está fora da altura do cone
    end

    r_ = sqrt(x_^2 + y_^2)  # Distância radial no plano
    r_cone = radius * (1 - (z_/height))  # Raio do cone no plano z

    return r_ > r_cone*1.05  # Retorna true se o ponto estiver fora do cone
end

function is_outside_sphere(x_::Float64, y_::Float64, z_::Float64, radius::Float64)
    r_ = sqrt(x_^2 + y_^2 + z_^2)  # Distância radial no plano
    return r_ > radius  # Retorna true se o ponto estiver fora do cone
end

function is_outside(type, a, b, c, radius, height=0)
    if type == "Cone"
        return is_outside_cone(a,b,c,radius,height)
    elseif type == "Sphere"
        return is_outside_sphere(a,b,c,radius)
    else
        return true
    end
end

norm(vec) = sqrt(sum(vec.^2))

struct charge
    r::Tuple{Float64, Float64, Float64}
    Q::Float64
end

function calc_E_field(s::surface, r::Tuple{Float64, Float64, Float64}, charges::Vector{charge}=charge[])
    E = [0.0, 0.0, 0.0]
    tol = 0.0
    A = area(s)
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

# function calc_potential(s::surface, r::Tuple{Float64, Float64, Float64})
#     V = 0.0
#     A = area(s)
#     for e in s.elements
#         rvec = [r[1] - e.center[1], r[2] - e.center[2], r[3] - e.center[3]]
#         q = e.area / A
#         ramp = norm(rvec)
#         if ramp > 0.0
#             V += q / ramp
#         end
#     end
#     return V
# end

function cone_surf_xyz(R, H, subdiv=10)
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
    return x, y, Z
end

R = 1.0
H = 2.0
surf = gen_cone(R, H, 20);
surf.Q = -2.0

println(length(surf.elements))
cargas = [
    charge((-1.5, 0.0, 1.0), 1.0),
    charge((1.5, 0.0, 1.0), 1.0),
    # charge((0.0, 0.0, -0.3), 1.0),
]
qsx = [e.r[1] for e in cargas]
qsy = [e.r[2] for e in cargas]
qsz = [e.r[3] for e in cargas]

x_range = range(-1.5, 1.5, length=8)
y_range = range(-1.5, 1.5, length=8)
z_range = range(-0.2, 2.5, length=8)
gr = [(x, y, z) for x in x_range for y in y_range for z in z_range]
vx, vy, vz, x, y, z, magnitudes = [], [], [], [], [], [], []

# surf = surface([], 0.0)

for (xi, yi, zi) in gr
    if is_outside("Cone", xi, yi, zi, R, H)
        E = calc_E_field(surf, (xi, yi, zi), cargas)
        push!(vx, E[1])
        push!(vy, E[2])
        push!(vz, E[3])
        push!(x, xi)
        push!(y, yi)
        push!(z, zi)
        push!(magnitudes, norm(E))
    end
end
magnitudes ./= maximum(magnitudes)

colormap = cm.viridis
colors = colormap(magnitudes)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection="3d")
if length(cargas) > 0
    ax.plot(qsx, qsy, qsz, "oC3", ms=5)
end
if surf.Q != 0.0
    ax.plot_surface(cone_surf_xyz(R, H, 30)..., alpha=1.0, color="black")
end
ax.quiver(x, y, z, vx, vy, vz, length=0.1, normalize=true, color=colors)
ax.set_aspect("equal")
plt.show()

# r0 = (1.1, 0.0, 0.0)
# v0 = calc_potential(surf, r0)

# function findV(s::surface, v0, r0, lista, d=0.1, maxit=10)
#     drange = [-1.0*abs(d), +1.0*abs(d)]
#     gri = [(r0[1]+dx, r0[2]+dy, r0[3]+dz) for dx in drange for dy in drange for dz in drange]
#     V = [calc_potential(s, p) for p in gri]
#     idx = argmin(abs.(V .- v0))
#     r = [gri[idx][1], gri[idx][2], gri[idx][3], V[idx]]
#     push!(r, V[idx])
#     push!(lista, r)
#     if length(lista) > maxit
#         return
#     else
#         return findV(s, v0, r, lista, d, maxit)
#     end
# end

# xyz = []
# findV(surf, v0, r0, xyz, 0.5, 500)
# x_ = [i[1] for i in xyz]
# y_ = [i[2] for i in xyz]
# z_ = [i[3] for i in xyz]
# V_ = [i[4] for i in xyz]
