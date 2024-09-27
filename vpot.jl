using DifferentialEquations, LinearAlgebra, Makie

# Estrutura para uma carga pontual
struct charge
    r::Vector{Float64}  # Posição
    Q::Float64          # Carga
end

# Função para calcular o campo elétrico em um ponto r
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


# Definir a direção perpendicular ao campo elétrico
function equipotential_direction!(du, u, p, t)
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

# Crie algumas cargas como exemplo
charges = [charge([1.0, 0.0, 0.0], 1.0), charge([-1.0, 0.0, 0.0], -1.0)]

dz = 0.2
x, y, z = [], [], []
zi = [0.0, ]

# tspan = (0.0, 10.0)  # Intervalo de "tempo" (parâmetro ao longo da superfície)
# st = [0.5, 0.5, 0.0]
# u0 = st  # Condição inicial (posição inicial)
# prob = ODEProblem(equipotential_direction!, u0, tspan, charges)

# # Resolver usando o solver padrão (Runge-Kutta 4)
# sol = solve(prob, Tsit5())  # Tsit5 é um método adaptativo de Runge-Kutta

# x = [r[1] for r in copy(sol.u)]
# y = [r[2] for r in copy(sol.u)]
# z = [r[3] for r in copy(sol.u)]

xc = [c.r[1] for c in charges]
yc = [c.r[2] for c in charges]
zc = [c.r[3] for c in charges]

function f(charges, xy, z, V0)
    return calc_V(charges, [xy[1], xy[2], z]) - V0
end

V0 = calc_V(charges, st0)

guess = [st0[1], st0[2]]

z = st0[3] + dz

result = nlsolve(xy -> f(charges, xy, z, V0), guess)

# Plotar a curva equipotencial resultante
st0 = [0.5, 0.5, 0.0]
if true
    fig = Figure()
    ax = Axis3(fig[1,1])
    for zi in range(-1.0, 1.0, 10)
        tspan = (0.0, 10.0)  # Intervalo de "tempo" (parâmetro ao longo da superfície)
        xi, yi = result.zero
        st = [0.5 - zi, 0.5, 0.0 + zi]

        V0 = calc_V(charges, st0)
        guess = [st0[1], st0[2]]
        result = nlsolve(xy -> f(charges, xy, zi, V0), guess)
        xi, yi = result.zero
        u0 = [xi, yi, zi]  # Condição inicial (posição inicial)
        prob = ODEProblem(equipotential_direction!, u0, tspan, charges)

        # Resolver usando o solver padrão (Runge-Kutta 4)
        sol = solve(prob, Tsit5())  # Tsit5 é um método adaptativo de Runge-Kutta

        x = [r[1] for r in copy(sol.u)]
        y = [r[2] for r in copy(sol.u)]
        z = [r[3] for r in copy(sol.u)]
        if zi == 0 || zi == 1.0
            println(calc_V(charges, sol.u[end]), "  $(sol.u[end])", "  ", calc_V(charges, sol.u[1]), "  $(sol.u[1])")
        end
        scatter!(ax, x, y, z, color=:blue, linewidth=2)
    end
    scatter!(ax, xc, yc, zc, color=:red)  # Posições das cargas
    display(fig)
end

using  NLsolve
