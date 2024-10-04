using DifferentialEquations, LinearAlgebra, GLMakie, NLsolve
using Meshes
using MeshIO, FileIO
using MeshBridge

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

function calc_E(charges::Vector{charge}, r::Vector{Float64}, tol=1e-9)
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

g_area(s::charged_surface) = length(s.elements) == 0 ? 0.0 : sum([e.area for e in s.elements])

function charged_surface_to_point_charges(s::charged_surface)
    A = sum([e.area for e in s.elements])
    Q = s.Q
    return [charge(e.centroid, Q * e.area/A) for e in s.elements]
end

R1 = 3
r1 = (-5.0, 0.0, 0.0)
esf1 = Meshes.Sphere(r1, R1)
msh = Meshes.discretize(esf1, Meshes.RegularDiscretization(30))
Q0 = 1.0
charged_surf = charged_surface(msh, Q0, 2)
surf_charges = charged_surface_to_point_charges(charged_surf)

R2 = 1.0
r2 = (5.0, 0.0, 0.0)
esf2 = Meshes.Sphere(r2, R2)
msh = Meshes.discretize(esf2, Meshes.RegularDiscretization(30))
charged_surf = charged_surface(msh, -Q0, 2)
surf_charges2 = charged_surface_to_point_charges(charged_surf)

# baloon = load("airballoon.obj")
# baloon = load("Ciclo2/airplane_model/11805_airplane_v2_L2.obj")
# baloon = load("12222_Cat_v1_l3.obj")

# mbaloon = convert(Meshes.Mesh, baloon)
# ba = charged_surface(mbaloon, -Q0, 2)
# qba = charged_surface_to_point_charges(ba)

# chs = [charge([1.0, 0.0, 0.0], 1.0), charge([-1.0, 0.0, 0.0], -1.0)]

charges = charge[]
append!(charges, surf_charges)
append!(charges, surf_charges2)
# append!(charges, qba)
# append!(charges, chs)

xc1 = [c.r[1] for c in charges[1:length(surf_charges)]]
yc1 = [c.r[2] for c in charges[1:length(surf_charges)]]
zc1 = [c.r[3] for c in charges[1:length(surf_charges)]]

xc = [c.r[1] for c in charges[length(surf_charges):end]]
yc = [c.r[2] for c in charges[length(surf_charges):end]]
zc = [c.r[3] for c in charges[length(surf_charges):end]]

function f(charges, xy, z, V0)
    return calc_V(charges, [xy[1], xy[2], z]) - V0
end

Base.GC.gc()

plot_flag = Bool(1)
calcE_flag = Bool(1)
calcV_flag = Bool(0)

if true
    if plot_flag
        fig = Figure()
        ax = Axis3(fig[1,1], aspect=(1,1,1))
    end
    # Linhas de Campo
    if calcE_flag
        dt = 2pi/15
        starts = [
            [r1[1]+ (R1+0.01)*cos(to), r1[2], r1[3] + (R1+0.01)*sin(to)] for to in range(0, 2pi-dt, step=dt) if ((abs(to-pi) > 0.001) && (to != 0.0))
        ]
        println(length(starts))
        append!(starts, [[r1[1] + (R1+0.01)*sin(to), r1[2] + (R1+0.01)*cos(to), r1[3]] for to in range(0, 2pi-dt, step=dt)])
        println(length(starts))
        Nd = 10
        ts = [
            Nd, Nd, Nd, Nd, Nd,
            Nd, Nd, Nd, Nd, Nd,
            Nd, Nd, Nd, Nd,
            Nd, Nd, Nd, Nd, Nd,
            Nd, Nd, Nd, Nd, Nd,
            Nd, Nd, Nd, Nd, Nd
        ]
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
                lines!(ax, x, y, z, color=:green, linewidth=2, alpha=1.0)
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

    if calcV_flag
        # Superficie Equipotencial
        st0 = [25.0, 0.0, 15.0]
        V0 = calc_V(charges, st0)
        guess = [st0[1], st0[2]]
        dz = 0.01
        z = st0[3] + dz
        result = nlsolve(xy -> f(charges, xy, z, V0), guess)
        Nff = 30
        for (ilo, zi) in enumerate(range(-25, 47, Nff))
            println("runing V line $(ilo)")
            tspan = (0.0, 150 * sin(pi*ilo/(Nff*1.2))^2)  # Intervalo de "tempo" (parâmetro ao longo da superfície)
            println(tspan)
            xi, yi = result.zero
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
                lines!(ax, x, y, z, color=:blue, linewidth=1, alpha=0.5)
            end
        end
        println("end V surface")
    end

    # scatter!(ax, xc, yc, zc, color=:red, markersize=5)  # Posições das cargas
    # scatter!(ax, xc1, yc1, zc1, color=:grey, markersize=5)  # Posições das cargas
    if plot_flag
        Meshes.viz!(esf1)
        Meshes.viz!(esf2)
        display(fig)
        d = 26
        xlims!(-d/2, d/2)
        ylims!(-d/2, d/2)
        zlims!(-d/2, d/2)
    end
end
