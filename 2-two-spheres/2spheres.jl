using DifferentialEquations, LinearAlgebra, GLMakie, NLsolve
using Meshes

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

R = 1.0/1.5
r0 = (-0.5, 0.5, 0.0)
esf = Meshes.Sphere(r0, R)
msh = Meshes.discretize(esf, Meshes.RegularDiscretization(20))
Q0 = 1.0
charged_surf = charged_surface(msh, Q0, 2)
surf_charges = charged_surface_to_point_charges(charged_surf)

R2 = 0.3/1.5
r02 = (1.0, 0.0, 0.0)
esf2 = Meshes.Sphere(r02, R2)
msh2 = Meshes.discretize(esf2, Meshes.RegularDiscretization(15))
charged_surf2 = charged_surface(msh2, -Q0/3, 2)
surf_charges2 = charged_surface_to_point_charges(charged_surf2)

R3 = 0.3/1.5
r03 = (0.0, -1.0, 0.0)
esf3 = Meshes.Sphere(r03, R2)
msh3 = Meshes.discretize(esf3, Meshes.RegularDiscretization(10))
charged_surf3 = charged_surface(msh3, -Q0/3, 2)
surf_charges3 = charged_surface_to_point_charges(charged_surf3)

charges = charge[]
append!(charges, surf_charges)
append!(charges, surf_charges2)
append!(charges, surf_charges3)

function f(charges, xy, z, V0)
    return calc_V(charges, [xy[1], xy[2], z]) - V0
end

fig = Figure(size=(800, 600))
ax = Axis3(fig[1,1], aspect=(1,1,1))
xlims!(-d/2, d/2)
ylims!(-d/2, d/2)
zlims!(-d/2, d/2)
Meshes.viz!(msh, alpha=1.0, color=:red)
Meshes.viz!(msh2, alpha=1.0, color=:blue)
Meshes.viz!(msh3, alpha=1.0, color=:blue)
save("scheme.png", fig)

Base.GC.gc()

plot_flag = true
calcE_flag = true
calcV_flag = true
quiver_flag = false
animation_flag = true

function Efunc(k)
    e_ = calc_E(charges, [k[1], k[2], k[3]])
    return GLMakie.Point3f(e_[1],e_[2],e_[3])
end

d = 5*R
if true
    if plot_flag
        fig = Figure(size=(800, 600))
        ax = Axis3(fig[1,1], aspect=(1,1,1))
        xlims!(-d/2, d/2)
        ylims!(-d/2, d/2)
        zlims!(-d/2, d/2)
        Meshes.viz!(msh, alpha=1.0, color=:red)
        Meshes.viz!(msh2, alpha=1.0, color=:blue)
        Meshes.viz!(msh3, alpha=1.0, color=:blue)
    end

    if calcV_flag
        # Superficie Equipotencial
        st0 = [0.35, -0.35, 0.0]
        V0 = calc_V(charges, st0)
        guess = [st0[1], st0[2]]
        z = st0[3]
        result = nlsolve(xy -> f(charges, xy, z, V0), guess)
        println(result)
        Nff = 50

        start_angle = pi/6
        n_frames = Nff
        ax.viewmode = :fit

        zz = range(-R*1.7, R*1.7, Nff)
        record(fig, "potential.mp4", 1:n_frames, framerate=15) do ilo
        # for (ilo, zi) in enumerate(range(-R*1.5, R*1.5, Nff))
            zi = zz[ilo]
            println("runing V line $(ilo)")
            tspan = (0.0, 20 * (0.3 + sin(pi*(ilo-1)/Nff))^3)  # Intervalo de "tempo" (parâmetro ao longo da superfície)
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
                save("potential/frame$(ilo).png", fig)
            end
        end
        println("end V surface")
    end

    # Linhas de Campo
    if calcE_flag
        dt = 2pi/8
        k = R+0.05
        starts = [
            [r0[1] + k*cos(to), r0[2] + k*sin(to), r0[3] + 0.0] for to in range(0, 2pi-dt, step=dt)
        ]
        append!(starts, [[r0[1] + k*cos(to), r0[2] + 0.0, r0[3] + k*sin(to)] for to in range(0, 2pi-dt, step=dt) if (to!=0) && (abs(to-pi)>1e-4)])
        k = (R+0.05)*cos(pi/4)
        h = (R+0.05)*sin(pi/4)
        append!(starts, [[r0[1] + k*cos(to), r0[2] + k*sin(to), r0[3] + h] for to in range(0, 2pi-dt, step=dt)])
        k = (R+0.05)*cos(-pi/4)
        h = (R+0.05)*sin(-pi/4)
        append!(starts, [[r0[1] + k*cos(to), r0[2] + k*sin(to), r0[3] + h] for to in range(0, 2pi-dt, step=dt)])
        # push!(starts, [r02[1]+R2+0.05, r02[2]-R2-0.05, r02[3]])
        # push!(starts, [r03[1]-R3-0.01, r03[2]+R3+0.01, r03[3]])
        Nd = 3.5
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
            name = "efield"
        else
            name = "full"
        end

        record(fig, name*".mp4", 1:n_frames, framerate=10) do frame
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
