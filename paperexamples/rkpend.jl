using DifferentialEquations

# m1=m2=l1=l2=1

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

# function f(u,p,t)
#
# end

u0 = [pi/2;0.;0;0]
tspan = (0.0,3600.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, RK4(),dt=0.01,adaptive=false)

myu=zeros(4,length(sol.u))
M = [zeros(2,2) for i=1:length(sol.u)]
T = zeros(length(sol.u))
V = zeros(length(sol.u))
E = zeros(length(sol.u))

for (i,vec) in enumerate(sol.u)
    for (j,el) in enumerate(vec)
        myu[j,i] = el
    end
    # M[i] = [2 cos(vec[1]-vec[2]);cos(vec[1]-vec[2]) 1]
    T[i] = .5*vec[3]^2 + .5*(vec[3]^2 + vec[4]^2 + 2*vec[3]*vec[4]*cos(vec[1]-vec[2]))#.5*vec[3:4]'*M[i]*vec[3:4]
    V[i] = -2*9.81*cos(vec[1]) - 9.81*cos(vec[2])#1*9.81*(.5*sin(vec[1]) + sin(vec[1])+.5*sin(vec[1]+vec[2]))
    E[i] = T[i]+V[i]
end

using Plots

# plot(sol)
