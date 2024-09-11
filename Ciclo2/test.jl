using PyCall

plt = pyimport("matplotlib.pyplot")
cm = pyimport("matplotlib.cm")

struct Surface
    points::Array{Tuple{Float64, Float64, Float64}, 1}
    x
    y
    z
    type::String
    function Surface(points, type="None")
        x, y, z = [], [], []
        for p in points
            push!(x, p[1])
            push!(y, p[2])
            push!(z, p[3])
        end
        return new(points, x, y, z, type)
    end
end

function gen_cone(radius, heigth, num=10, num2=100)
    num = Int(ceil(num))
    num2 = Int(ceil(num2))
    vertices = []
    for i in 0:num-1
        r = radius * (1 - i / (num - 1))
        z = i * heigth / (num - 1)
        num_ = Int(ceil(num2 * r / radius) + 1)
        if num_ != 1
            rd = rand()*2pi
            theta = range(rd, rd+2pi, num_)
            x = r * sin.(theta)
            y = r * cos.(theta)
            for j in 1:num_
                push!(vertices, (x[j],y[j],z))
            end
        else
            push!(vertices, (0.0, 0.0, heigth))
        end
    end
    Abase = pi * radius^2
    Alat = pi * radius * sqrt(radius^2 + heigth^2)
    leng1 = length(vertices) * 1.0
    # println(leng1)
    siz = leng1 * Abase / Alat
    # println(siz)
    siz = Int(floor(sqrt(siz*1.2)))
    for x in range(-radius, radius, siz)
        for y in range(-radius, radius, siz)
            if x^2 + y^2 < (radius*(1-(1/num)))^2
                push!(vertices, (x,y,0.0))
            end
        end
    end
    # println(length(vertices)-leng1)
    return Surface(vertices, "Cone")
end

function gen_sphere(radius, num=100)
    num = Int(ceil(num))
    vertices = []
    for phi in range(-pi, pi, num)
        z = radius * sin(phi)
        r = sqrt(radius^2 - z^2)
        num_ = Int(floor(num * r / radius) + 1)
        rd = rand()*2pi
        for theta in range(rd, rd+2pi, num_)
            x = r * sin(theta)
            y = r * cos(theta)
            push!(vertices, (x,y,z))
        end
    end
    return Surface(vertices, "Sphere")
end

function is_outside_cone(x_::Float64, y_::Float64, z_::Float64, radius::Float64, height::Float64)
    println("z = ", z_)
    if (z_ < 0.0) || (z_ > height)
        return true  # Ponto está fora da altura do cone
    end

    r_ = sqrt(x_^2 + y_^2)  # Distância radial no plano
    r_cone = radius * (1 - (z_/height))  # Raio do cone no plano z

    return r_ > r_cone*1.1  # Retorna true se o ponto estiver fora do cone
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

function calc_E_field(points::Array{Tuple{Float64, Float64, Float64}, 1}, r::Tuple{Float64, Float64, Float64})
    Ex, Ey, Ez = 0.0, 0.0, 0.0
    for p in points
        rvec = (r[1]-p[1], r[2]-p[2], r[3]-p[3])
        ramp = sqrt(rvec[1]^2 + rvec[2]^2 + rvec[2]^2)
        # if ramp != 0.0
        Ex += (rvec[1] / (ramp^3))
        Ey += (rvec[2] / (ramp^3))
        Ez += (rvec[3] / (ramp^3))
        # end
    end
    return (Ex, Ey, Ez)
end

R, H = 1.0, 2.0
N = 200
N2 = N
# surf = gen_cone(R, H, N, N2);
surf = gen_sphere(R, N)

# Generate a grid of points (x, y, z)
x_range = range(-3, 3, length=25)
y_range = range(-3, 3, length=25)
# z_range = range(-1.1, 1.1, length=7)
z_range = [2.0,]

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
    if is_outside(surf.type, xi, yi, zi, R, H)
        vxi, vyi, vzi = calc_E_field(surf.points, (xi, yi, zi))
        push!(vx, vxi)
        push!(vy, vyi)
        push!(vz, vzi)
        push!(x, xi)
        push!(y, yi)
        push!(z, zi)
        push!(magnitudes, sqrt(vxi^2 + vyi^2 + vzi^2))
    end
end

magnitudes ./= maximum(magnitudes)

plt.plot(magnitudes)
plt.show()

colormap = cm.viridis
colors = colormap(magnitudes)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection="3d")
# ax.plot_trisurf(surf.x, surf.y, surf.z)
ax.plot(surf.x, surf.y, surf.z, ".")
ax.quiver(x, y, z, vx, vy, vz, length=0.2, normalize=true, color=colors)
ax.set_aspect("equal")
# ax.set_zlim(0, 0.02)
plt.show()


# # # for xi in x_range
# # #     for yi in y_range
# # #         for zi in z_range
# # #             c = "C1"
# # #             alpha=1.0
# # #             if is_outside(surf.type, xi, yi, zi, R, H)
# # #                 c = "C0"
# # #                 alpha = 0.2
# # #             end
# # #             ax.plot([xi], [yi], [zi], ".", color=c, alpha=alpha)
# # #         end
# # #     end
# # # end
