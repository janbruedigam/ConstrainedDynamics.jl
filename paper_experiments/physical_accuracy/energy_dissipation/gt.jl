using DifferentialEquations

# mass1 = mass2 = link_length1 = link_length2 = 1
J = 0.3333

# Pendulum differential equations
function f(u,p,t)
    [
    u[2]
    -9.81*3/2*sin(u[1]) - u[2]*3/2
    ]
end

u0 = [π/2;0]
timespan = (0.0,10.0)
problem = ODEProblem(f,u0,timespan)
sol = solve(problem, Tsit5(),dt=0.0001,adaptive=false)

energy = -9.81*cos.(getindex.(sol.u,1))/2 + 0.5*J*getindex.(sol.u,2).^2

# plot!(sol.t,getindex.(sol.u,1))