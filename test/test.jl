using ConstrainedDynamics
using LinearAlgebra
using Rotations

length1 = 1.0
width, depth = 1.0, 1.0
box = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
box2 = Box(width*2, depth*2, length1*2, length1*2, color = RGBA(1., 1., 0.))

# Links
origin = Origin{Float64}()
link1 = Body(box)
link2 = Body(box2)

# Constraints
joint1 = EqualityConstraint(Spherical(origin, link1))
joint2 = EqualityConstraint(Revolute(origin, link2, [1.;0;0]))

links = [link1;link2]
constraints = [joint1;joint2]

mech = Mechanism(origin, links, constraints)

xd=zeros(6)
vd=zeros(6)
Fd=zeros(6)
qd=Quaternion{Float64}()
qd=[qd;qd]
ωd=zeros(6)
ωd[1] = 1.
τd=zeros(6)

A, B, G = linearizeMechanism(mech, xd, vd, Fd, qd, ωd, τd)