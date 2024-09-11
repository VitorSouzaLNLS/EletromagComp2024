using PyCall

plt = pyimport("matplotlib.pyplot")
cm = pyimport("matplotlib.cm")
np = pyimport("numpy")

struct element
    center
    area
end

struct surface
    elements::Vector{element}
    surface(elems) = new(elems)
end

area(s::surface) = sum([e.area for e in s.elements])

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
    if (z_ < 0.0) || (z_ > height)
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

function calc_E_field(s::surface, r::Tuple{Float64, Float64, Float64})
    E = [0.0, 0.0, 0.0]
    A = area(s)
    for e in s.elements
        rvec = [r[1] - e.center[1], r[2] - e.center[2], r[3] - e.center[3]]
        q = e.area / A
        ramp = norm(rvec)^3
        if ramp > 0.0
            E = E .+ (q .* rvec ./ ramp)
        end
    end
    return E
end

R = 1.0
H = 2.0
surf = gen_cone(R, H, 10);
println(length(surf.elements))

sx = [e.center[1] for e in surf.elements]
sy = [e.center[2] for e in surf.elements]
sz = [e.center[3] for e in surf.elements]

x_range = range(-1.3, 1.3, length=7)
y_range = range(-1.3, 1.3, length=7)
z_range = range(-0.2, 2.2, length=7)
# x_range = [0.0,]
# y_range = [0.0,]
# z_range = [1.2,]

gr = [(x, y, z) for x in x_range for y in y_range for z in z_range]

# Calculate vectors for each point
vx = Float64[]
vy = Float64[]
vz = Float64[]
x = Float64[]
y = Float64[]
z = Float64[]
magnitudes = Float64[]

for (xi, yi, zi) in gr
    if is_outside("Cone", xi, yi, zi, R, H)
        E = calc_E_field(surf, (xi, yi, zi))
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

plt.plot(magnitudes)
plt.show()

colormap = cm.viridis
colors = colormap(magnitudes)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection="3d")
# ax.plot_trisurf(sx, sy, sz)
ax.plot(sx, sy, sz, ".C1", alpha=0.5)
# ax.quiver(x, y, z, vx, vy, vz, length=0.1, normalize=false, color="black")
ax.quiver(x, y, z, vx, vy, vz, length=0.1, normalize=true, color=colors)
ax.set_aspect("equal")
# ax.set_zlim(0, 0.02)
plt.show()
