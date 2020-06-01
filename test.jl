using ConstrainedDynamics
using LinearAlgebra
using Rotations

length1 = 1.0
width, depth = 1.0, 1.0
box = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
box2 = Box(width*2, depth*2, length1*2, length1*2, color = RGBA(1., 1., 0.))

# Links
origina = Origin{Float64}()
originb = Origin{Float64}()
link1a = Body(box)
link2a = Body(box2)
link1b = Body(box)
link2b = Body(box2)

# Constraints
joint1 = EqualityConstraint(Spherical(originb, link1b))
joint2 = EqualityConstraint(Revolute(originb, link2b, [1.;0;0]))

linksa = [link1a;link2a]
linksb = [link1b;link2b]
constraints = [joint1;joint2]
shapes = [box]


mecha = Mechanism(origina, linksa, shapes = shapes)

mechb = Mechanism(originb, linksb, constraints, shapes = shapes)

# Q = diagm(ones(12))
# Q[2,2] = 2
# Q[3,3] = 0.1
# Q[7,7] = 2
# R = diagm(ones(6))
# R[1] = 100


xd=zeros(6)
vd=zeros(6)
Fd=zeros(6)
qd=Quaternion{Float64}()
qd=[qd;qd]
ωd=zeros(6)
ωd[1] = 1.
τd=zeros(6)
    
A, B = ConstrainedDynamics.linearizeSystem(mecha, xd, vd, Fd, qd, ωd, τd)
G = ConstrainedDynamics.linearizeConstraints(mechb, xd, vd, Fd, qd, ωd, τd)