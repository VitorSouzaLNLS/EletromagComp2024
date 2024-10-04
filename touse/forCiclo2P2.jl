using DifferentialEquations


struct charge
    Q::Float64
    r::Vector{Float64}
    charge(q::Number=0.0, r) = new(Float64(q), Float64.([r[1], r[2], r[3]]))
end



struct surface_fragment
    centroid::Vector{Float64}
    area::Float64
end



mutable struct charged_surface
    Q::Float64=0.0
    elements::Vector{surface_fragment}
end



function charged_surface(m::Mesh, Q::Number)
    elems = surface_fragment[]
    for fragment in faces(m, 3)
        push!(elems, surface_fragment(coordinates(centroid(fragment)), area(fragment)))
    end
    return surface(elems, Q)
end


area(s::charged_surface) = length(s.elements) == 0 ? 0.0 : sum([e.area for e in s.elements])


function charged_surface_to_point_charges(s::surface)
    A = sum([e.area for e in s.elements])
    Q = s.Q
    return Vector{charge}([charge(Q * e.area/A, e.centroid) for e in s.elements])
end



function calc_E(charges::Vector{charge}=charge[], r::Vector{Float64})
    E = [0.0, 0.0, 0.0]
    for e in charges
        rvec = r .- e.r
        ramp = norm(rvec)
        if ramp > tol
            E .+= (e.Q .* rvec ./ ramp^3)
        end
    end
    return E
end



function calc_V(charges::Vector{charge}=charge[], r::Vector{Float64})
    V = 0.0
    for e in charges
        rvec = r .- e.r
        ramp = norm(rvec)
        if ramp > tol
            V += (e.Q / ramp)
        end
    end
    return V
end



function Efield_line_EDO!(dr, r, charges::Vector{charge}=charge[], s)
    E = calc_E(charges, r)
    dr[:] = E ./ norm(E)
end



function get_Efield_line(charges::Vector{charge}=charge[], r0::Vector{Float64})
    s_range = (0.0, 30.0)
    edo = ODEProblem(Efield_line_EDO!, r0, s_range, charges)
    return solve(edo)
end


function Vfield_line_EDO!(dr, r, charges::Vector{charge}=charge[], s)
    E = calc_E(charges, r)
    v = [0.0, 0.0, 1.0]
    n = cross(E, v)
    dr[:] = n ./ norm(n)
end



function get_Vfield_line(charges::Vector{charge}=charge[], r0::Vector{Float64})
    s_range = (0.0, 30.0)
    edo = ODEProblem(Vfield_line_EDO!, r0, s_range, charges)
    return solve(edo)
end




function func(charges, x_y, z, V0)
    return calc_V(charges, [x_y[1], x_y[2], z]) - V0
end

V0 = calc_V(charges, [x0, y0, z0])
guess = [x0, y0]
result = nlsolve(x_y -> func(charges, x_y, z0 + dz, V0), guess)
xi, yi = result.zero
r = [xi, yi, z0 + dz]
