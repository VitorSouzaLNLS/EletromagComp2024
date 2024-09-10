using PyCall

plt = pyimport("matplotlib.pyplot")


mutable struct Circuit
    R::Float64
    L::Float64
    C::Float64
    V0::Float64
    freq::Float64
    function Circuit(r=1, l=1, c=1, v0=1, f=0)
        return new(r, l, c, v0, f)
    end
end

#Funcao para definir o sistema de EDOs
function rlc(u, circuit, ti)
    I, Vc = u
    V = circuit.V0
    if circuit.freq !== 0.0
        V *= sinpi(2*circuit.freq * ti)
    end
    return [
        (1/circuit.L) * (V - circuit.R*I - (1/circuit.C)*Vc),  #dI/dt
        I  #dVc/dt
    ]
end

#------------------------------------------------------

#Funcao para resolver o sistema de EDOs usando RK4
function rk4(f::Function, u0::Vector{Float64}, tspan::Tuple{Float64, Float64}, params::Circuit, dt::Float64=0.01)
    t0, tf = tspan
    N = Int((tf - t0) / dt) + 1
    t = t0:dt:tf
    u = zeros(N, length(u0))
    u[1, :] = u0

    for i in 1:(N-1)
        ti = t[i]
        k1 = dt * f(u[i, :], params, ti)
        k2 = dt * f(u[i, :] .+ 0.5 * k1, params, ti)
        k3 = dt * f(u[i, :] .+ 0.5 * k2, params, ti)
        k4 = dt * f(u[i, :] .+ k3, params, ti)
        u[i+1, :] = u[i, :] .+ (k1 .+ 2*k2 .+ 2*k3 .+ k4) ./ 6
    end

    return t, u
end

#------------------------------------------------------

#Definicao dos parametros

R, L, C = 2, 0.5, 0.2
V0 = 10

malha = Circuit(R, L, C, V0)

#Condicoes iniciais
u0 = [0.0, 0.0]  #[I(0), Vc(0)]
tspan = (0.0, 20.0)  #Intervalo de tempo
dt = 1e-4  #Passo de tempo

#------------------------------------------------------

#Execucao o metodo RK4
Res = []
malha.freq = 0
Rs = [0]
labels = ["Crítico"]
labels = [""]
# Rs = [-1, 0, +1]
# labels = ["Sub Amortecido", "Crítico", "Super Amortecido"]

for i in Rs
    malha.R = R + i
    println(malha.R)
    t, u = rk4(rlc, u0, tspan, malha, dt)
    push!(Res, (t, u))
end

ResAC = []
malha.freq = 1/sqrt(L*C) / (2pi)
for i in Rs
    malha.R = R + i
    println(malha.R)
    t, u = rk4(rlc, u0, tspan, malha, dt)
    push!(ResAC, (t, u))
end

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
# for (i,r) in enumerate(Res)
#     ax[1].plot(r[1], r[2][:,1], label=labels[i])
#     ax[2].plot(r[1], r[2][:,2], label=labels[i])

# end
for (i,r) in enumerate(ResAC)
    ax[1].plot(r[1], r[2][:,1], label=labels[i])
    ax[2].plot(r[1], r[2][:,2], label=labels[i])
end
ax[1].set_xlabel("Tempo [s]")
ax[1].set_ylabel("Corrente [A]")
ax[1].legend()
ax[2].set_xlabel("Tempo [s]")
ax[2].set_ylabel("Tensão no Capacitor [V]")
ax[2].legend()

fig.tight_layout()
plt.show()
