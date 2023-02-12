using ConstrainedDynamics
using ConstrainedDynamicsVis
using BenchmarkTools

# Parameters
ex = [1.;0.;0.]

h = 1.
r = .05

vert11 = [0.;0.;h / 2]
vert12 = -vert11

# Initial orientation
phi = pi / 4
q1 = Quaternion(RotX(phi))

N = 7

origin = Origin{Float64}()
links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:N]

# Constraints
jointb1 = EqualityConstraint(Floating(origin, links[1]))
if N>1
    constraints = [jointb1;[EqualityConstraint(Spherical(links[i - 1], links[i]; p1=vert12, p2=vert11)) for i = 2:N]]
else
    constraints = [jointb1]
end

fricsandineqs1 = Friction(links[1], [0;0;1.0], 0.4; p = vert11)
frics1 = fricsandineqs1[1]
ineqcs1 = fricsandineqs1[2]
fricsandineqs2 = [Friction(links[i], [0;0;1.0], 0.4; p = vert12) for i=1:N]
frics2 = getindex.(fricsandineqs2,1)
ineqcs2 = vcat(getindex.(fricsandineqs2,2)...)

mech = Mechanism(origin, links, constraints, [ineqcs1;ineqcs2], [frics1;frics2])
# setPosition!(origin,links[1],p2 = [0;0.0;0.5],Δq = Quaternion(RotX(pi/2)*RotY(pi/4)))
setPosition!(origin,links[1],p2 = [0;-2.0;0.5],Δq = Quaternion(RotX(pi/2)*RotZ(pi/4)))
setVelocity!(links[1],v=randn(3)*2,ω=randn(3)*2)
for j=2:N
    if iseven(j)
        setPosition!(links[j-1],links[j],p1 = vert12, p2 = vert11,Δq = Quaternion(RotY(-pi/2)))
    else
        setPosition!(links[j-1],links[j],p1 = vert12, p2 = vert11,Δq = Quaternion(RotY(pi/2)))
    end  
    setVelocity!(links[j],v=randn(3)*2,ω=randn(3)*2)
end

steps = Base.OneTo(1000)
storagesnake = Storage(steps,length(mech.bodies))

storagesnake = simulate!(mech, storagesnake, record = true, ε=1e-6)
visualize(mech, storagesnake)