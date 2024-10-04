using DifferentialEquations, LinearAlgebra, GLMakie, NLsolve
using Meshes
using MeshBridge
using MeshIO, FileIO

struct charge
    r::Vector{Float64}  # Posição
    Q::Float64          # Carga
end

struct surface_fragment
    centroid::Vector{Float64}
    area::Float64
end

mutable struct charged_surface
    Q::Float64
    elements::Vector{surface_fragment}
end

function calc_E(charges::Vector{charge}, r::Vector, tol=1e-9)
    E = zeros(3)  # Vetor do campo elétrico (inicialmente zero)
    for e in charges
        rvec = r .- e.r
        ramp = norm(rvec)
        if ramp > tol
            E += e.Q * rvec / ramp^3  # Fórmula do campo elétrico
        end
    end
    return E
end

function calc_V(charges::Vector{charge}, r::Vector{Float64}, tol=1e-9)
    V = 0.0 # Vetor do campo elétrico (inicialmente zero)
    for e in charges
        rvec = r .- e.r
        ramp = norm(rvec)
        if ramp > tol
            V += e.Q  / ramp  # Fórmula do campo elétrico
        end
    end
    return V
end

function E_direction!(du, u, p, t)
    charges = p
    E = calc_E(charges, u)
    if norm(E) > 1e-9
        # Produto vetorial com vetor arbitrário para obter direção ortogonal
        du[:] = E / norm(E)
    end
end

function V_direction!(du, u, p, t)
    charges = p
    E = calc_E(charges, u)
    if norm(E) > 1e-9
        # Produto vetorial com vetor arbitrário para obter direção ortogonal
        v = [0.0, 0.0, 1.0]
        du[:] = cross(E, v) / norm(cross(E, v))
    else
        v = [0.0, 0.0, 1.0]
        du[:] = v  # Direção arbitrária se o campo for zero
    end
end

function charged_surface(m::Meshes.Mesh, Q::Number, discretizationtype=2)
    elems = surface_fragment[]
    for fragment in faces(m, discretizationtype)
        push!(elems, surface_fragment(coordinates(centroid(fragment)), area(fragment)))
    end
    return charged_surface(Q, elems)
end

function charged_surface_to_point_charges(s::charged_surface)
    A = sum([e.area for e in s.elements])
    Q = s.Q
    return [charge(e.centroid, Q * e.area/A) for e in s.elements]
end

Q0 = 1.0
R = 1.5
r1 = (-20.0, 15.0, 40.0)
esf1 = Meshes.Sphere(r1, R)
r2 = (25.0, 0.0, 40.0)
esf2 = Meshes.Sphere(r2, R)
r3 = (0.0, -25.0, 60.0)
esf3 = Meshes.Sphere(r3, R)

catload = load("/home/vitor-d/usp/EletromagComp2024/touse/12222_Cat_v1_l3.obj")
catmesh = convert(Meshes.Mesh, catload)
catsurf = charged_surface(catmesh, Q0, 2)
catcharges = charged_surface_to_point_charges(catsurf)

charges = charge[]
# append!(charges, surf_charges)
append!(charges, catcharges)
append!(charges, charge[charge([r1[1], r1[2], r1[3]], -Q0), charge([r2[1], r2[2], r2[3]], -Q0), charge([r3[1], r3[2], r3[3]], -Q0)])


function f(charges, xy, z, V0)
    return calc_V(charges, [xy[1], xy[2], z]) - V0
end

d = 100
ddd = 25

fig = Figure(size=(800, 600))
ax = Axis3(fig[1,1], aspect=(1,1,1))
xlims!(-d/3, d/3)
ylims!(-d/3, d/3)
zlims!(-d/3 + ddd, d/3 + ddd)
Meshes.viz!(catmesh, alpha=1.0, color=:cyan)
Meshes.viz!(esf1, alpha=1.0, color=:gray)
Meshes.viz!(esf2, alpha=1.0, color=:gray)
Meshes.viz!(esf3, alpha=1.0, color=:gray)
save("4-cat-and-sphere/1-mesh.png", fig)

# fig = Figure(size=(800, 600))
# ax = Axis3(fig[1,1], aspect=(1,1,3/2))
# xlims!(-d/2, d/2)
# ylims!(-d/2, d/2)
# zlims!(-d/2 + ddd, d/2 + ddd)
# Meshes.viz!(catmesh, alpha=1.0, color=:gray)
# scatter!(ax, xc[1:3:end], yc[1:3:end], zc[1:3:end], color=:red, markersize=5)
# save("2-mesh-and-charges.png", fig)

# fig = Figure(size=(800, 600))
# ax = Axis3(fig[1,1], aspect=(1,1,3/2))
# xlims!(-d/2, d/2)
# ylims!(-d/2, d/2)
# zlims!(-d/2 + ddd, d/2 + ddd)
# # Meshes.viz!(catmesh, alpha=1.0, color=:gray)
# scatter!(ax, xc, yc, zc, color=:red, markersize=5)
# save("3-charges.png", fig)

Base.GC.gc()

meansp = [(r1[1]+r2[1]+r3[1])/3 - 4, (r1[2]+r2[2]+r3[2])/3 , (r1[3]+r2[3]+r3[3])/3 ]

run_flag = false
plot_flag = true
calcE_flag = true
calcV_flag = false
quiver_flag = false
animation_flag = true

function Efunc(k)
    e_ = calc_E(charges, [k[1], k[2], k[3]])
    return GLMakie.Point3f(e_[1],e_[2],e_[3])
end

# println([i for i in 1:3_000:length(catcharges)])

k = [1]
for i in eachindex(catcharges[1:length(catcharges)])
    if all([(norm(catcharges[i].r .- catcharges[idx].r) > 6.0) for idx in k])
        push!(k, i)
    end
end
println(length(k))
selected = k

d = 100
ddd = 25
if run_flag
    if plot_flag
        fig = Figure(size=(800, 600))
        ax = Axis3(fig[1,1], aspect=(1,1,1))
        xlims!(-d/3, d/3)
        ylims!(-d/3, d/3)
        zlims!(-d/3 + ddd, d/3 + ddd)
        Meshes.viz!(catmesh, alpha=1.0, color=:cyan)
        Meshes.viz!(esf1, alpha=1.0, color=:gray)
        Meshes.viz!(esf2, alpha=1.0, color=:gray)
        Meshes.viz!(esf3, alpha=1.0, color=:gray)
    end

    if calcV_flag
        # Superficie Equipotencial
        st0 = meansp
        V0 = calc_V(charges, st0)
        guess = [st0[1], st0[2]]
        z = st0[3]
        result = nlsolve(xy -> f(charges, xy, z, V0), guess)
        println(result)
        Nff = 50

        start_angle = pi/6
        n_frames = Nff
        ax.viewmode = :fit

        # zz = range(-ddd+ddd/2.8, 2*ddd + ddd/18, Nff)
        zz = range(r1[3]-13, r3[3]+13, Nff)
        record(fig, "4-cat-and-sphere/potential.mp4", 1:n_frames, framerate=15) do ilo
        # for (ilo, zi) in enumerate(range(-R*1.5, R*1.5, Nff))
            zi = zz[ilo]
            println("runing V line $(ilo)")
            if zi <= meansp[3]
                tspan = (0.0, 100.0)
            else
                tspan = (0.0, 100.0)
            end
            println(tspan)
            # xi, yi = nlsolve(xy -> f(charges, xy, z, V0), guess).zero
            st = [st0[1] - zi, st0[2], zi]

            println("calculating x y")
            V0 = calc_V(charges, st0)
            guess = [st0[1], st0[2]]
            result = nlsolve(xy -> f(charges, xy, zi, V0), guess)
            xi, yi = result.zero
            u0 = [xi, yi, zi]  # Condição inicial (posição inicial)
            println("calculated")
            prob = ODEProblem(V_direction!, u0, tspan, charges)

            # Resolver usando o solver padrão (Runge-Kutta 4)
            sol = solve(prob, BS3())  # Tsit5 é um método adaptativo de Runge-Kutta

            x = [r[1] for r in copy(sol.u)]
            y = [r[2] for r in copy(sol.u)]
            z = [r[3] for r in copy(sol.u)]
            # if zi == 0 || zi == 1.0
            #     println(calc_V(charges, sol.u[end]), "  $(sol.u[end])", "  ", calc_V(charges, sol.u[1]), "  $(sol.u[1])")
            # end
            # scatter!(ax, x, y, z, color=:blue, linewidth=2)
            println("done solving V: line = $(ilo)")
            if plot_flag
                lines!(ax, x, y, z, color=:orange, linewidth=1, alpha=0.5)
                # ax.azimuth[] = start_angle + 2pi * frame / n_frames
                save("4-cat-and-sphere/potential/frame$(ilo).png", fig)
            end
        end
        println("end V surface")
    end

    # Linhas de Campo
    if calcE_flag
        starts = [catcharges[s].r for s in selected]
        Nd = 22.0
        ts = Float64[Nd for _ in eachindex(starts)]
        pvt = 0
        for (ij, st) in enumerate(starts)
            println("runing E line $(ij) --- $(ts[ij])")
            tspan = (0.0, ts[ij])  # Intervalo de "tempo" (parâmetro ao longo da superfície)

            u0 = [st[1], st[2], st[3]]  # Condição inicial (posição inicial)
            prob = ODEProblem(E_direction!, u0, tspan, charges)

            # Resolver usando o solver padrão (Runge-Kutta 4)
            sol = solve(prob, BS3())
            # sol = solve(prob, Tsit5())

            x = [r[1] for r in copy(sol.u)]
            y = [r[2] for r in copy(sol.u)]
            z = [r[3] for r in copy(sol.u)]

            # scatter!(ax, x, y, z, color=:blue, linewidth=2)
            if plot_flag
                lines!(ax, x, y, z, color=:black, linewidth=2, alpha=1.0)
            end
            # pvt += 1
            # if pvt == 5
            #     pvt = 0
            #     ts += 5
            # end
            println("E line $(ij) done")
        end
        println("end E lines")
    end


    if quiver_flag
        spec_ = 10
        xs = [xi for xi in range(-d/3, d/3, spec_)]
        ys = [yi for yi in range(-d/2, d/3, spec_)]
        zs = [zi for zi in range(-d*0.1, d*1.1, spec_)]
        grid = [[xi, yi, zi] for xi in range(-d/2, d/2, spec_) for yi in range(-d/2, d/2, spec_) for zi in range(-d/3, 2*d/3, spec_)]

        # mags = []
        # for p in grid
        #     e = calc_E(charges, p)
        #     # push!(Et, e)
        #     # push!(ex, e[1])
        #     # push!(ey, e[1])
        #     # push!(ez, e[3])
        #     push!(mags, norm(e))
        #     # push!(rx, p[1])
        #     # push!(ry, p[1])
        #     # push!(rz, p[3])
        # end
        tryfunc(x) = GLMakie.Point3f(x[1], x[2], x[3])

        println(length(xs))
        if plot_flag
            arrows!(ax, xs, ys, zs, Efunc, normalize=true, arrowsize=0.5, lengthscale=3, linewidth=0.3, alpha=0.8)
            # arrows!(ax, grid, Et)
        end
    end

    # if plot_flag
    #     Meshes.viz!(msh, alpha=1.0, color=:red)
    #     # scatter!(ax, xc, yc, zc, color=:red, markersize=5)  # Posições das cargas
    #     # lines!(ax, xc, yc, zc, color=:black)  # Posições das cargas

    #     xlims!(-d/2, d/2)
    #     ylims!(-d/2, d/2)
    #     zlims!(-d/2, d/2)
    # end

    if animation_flag
        start_angle = pi/6
        n_frames = 60
        ax.viewmode = :fit

        if !calcV_flag
            name = "4-cat-and-sphere/efield"
        else
            name = "4-cat-and-sphere/full"
        end

        record(fig, name*".mp4", 1:n_frames, framerate=15) do frame
            println(frame)
            ax.azimuth[] = start_angle + 2pi * frame / n_frames
        end

        for frame in range(1, n_frames)
            ax.azimuth[] = start_angle + 2pi * frame / n_frames
            save(name*"/frame$(frame).png", fig)
        end
    end

    display(fig)
end
