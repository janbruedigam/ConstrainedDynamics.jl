using DifferentialEquations

# mass1 = mass2 = link_length1 = link_length2 = 1

# Pendulum differential equations
function f(u,p,t)
    α1 = 0.5*cos(u[1]-u[2])
    α2 = cos(u[1]-u[2])

    f1 = -0.5*u[4]^2*sin(u[1]-u[2])-9.81*sin(u[1])
    f2 = u[3]^2*sin(u[1]-u[2])-9.81*sin(u[2])

    g1 = (f1-α1*f2)/(1-α1*α2)
    g2 = (-α2*f1+f2)/(1-α1*α2)

    [
    u[3]
    u[4]
    g1
    g2
    ]
end

u0 = [π/2;π/2;0;0]
timespan = (0.0,3600.0)
problem = ODEProblem(f,u0,timespan)
sol = solve(problem, Heun(),dt=0.01,adaptive=false)

T = zeros(length(sol.u)) # Kinetic energy
V = zeros(length(sol.u)) # Potential energy
energy = zeros(length(sol.u))

for (i,vec) in enumerate(sol.u)
    T[i] = .5*vec[3]^2 + .5*(vec[3]^2 + vec[4]^2 + 2*vec[3]*vec[4]*cos(vec[1]-vec[2]))
    V[i] = -2*9.81*cos(vec[1]) - 9.81*cos(vec[2])
    energy[i] = T[i]+V[i]
end