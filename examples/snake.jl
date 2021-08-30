using ConstrainedDynamics
using ConstrainedDynamicsVis
CD = ConstrainedDynamics

# Parameters
ex = [1.;0.;0.]

h = 1.
r = .05

vert11 = [0.;0.;h / 2]
vert12 = -vert11

# Initial orientation
phi = pi / 4
q1 = UnitQuaternion(RotX(phi))

# Links
N = 5

origin = Origin{Float64}()
links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:N]

# Constraints
jointb1 = EqualityConstraint(Floating(origin, links[1]))
if N>1
    constraints = [jointb1;[EqualityConstraint(Revolute(links[i - 1], links[i], ex; p1=vert12, p2=vert11)) for i = 2:N]]
else
    constraints = [jointb1]
end

fricsandineqs1 = [Friction(links[i], [0;0;1.0], 0.2; p = vert11) for i=1:N]
frics1 = getindex.(fricsandineqs1,1)
ineqcs1 = vcat(getindex.(fricsandineqs1,2)...)
fricsandineqs2 = [Friction(links[i], [0;0;1.0], 0.2; p = vert12) for i=1:N]
frics2 = getindex.(fricsandineqs2,1)
ineqcs2 = vcat(getindex.(fricsandineqs2,2)...)

mech = Mechanism(origin, links, constraints, [ineqcs1;ineqcs2], [frics1;frics2])
setPosition!(origin,links[1],p2 = [0;0.0;0.5],Δq = UnitQuaternion(RotX(pi/2)))
previd = links[1].id
for body in Iterators.drop(mech.bodies, 1)
    global previd
    setPosition!(ConstrainedDynamics.getbody(mech, previd), body, p1 = vert12, p2 = vert11)
    previd = body.id
end

storage = simulate!(mech, 10., record = true)
visualize(mech, storage)

# storage = Storage{Float64}(steps,length(mech.bodies))
# @btime simulate!($mech, $steps, $storage)