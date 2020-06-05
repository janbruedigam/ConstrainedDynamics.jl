using ConstrainedDynamics
using StaticArrays

# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 1.0
x, z = .1, .1
b1 = Box(x, l1, l1, l1 * 1, color = RGBA(1., 1., 0.))
b2 = Box(x, l1, z, l1 * 1, color = RGBA(1., 0., 0.))

vert11 = [0;l1 / 2;l1 / 2]
vert12 = -vert11

vert21 = [0;l1 / 2;0.]
vert22 = -vert21


# Links
origin = Origin{Float64}()
link1 = Body(b1)
link2 = Body(b2)

# Constraints
# joint0to1 = EqualityConstraint(OriginConnection(origin,link1))
# joint0to1 = EqualityConstraint(Fixed(origin,link1))
joint0to1 = EqualityConstraint(Revolute(origin, link1, ex))
joint1to2 = EqualityConstraint(Cylindrical(link1, link2, -ey; p1=vert12, p2=vert21))


links = [link1; link2]
constraints = [joint0to1;joint1to2]
shapes = [b1;b2]

function joint_force_control!(mechanism, k)
    F = SVector{2,Float64}(0.1,0.)
    setForce!(mechanism, joint1to2, F)
    return
end

mech = Mechanism(origin, links, constraints, shapes = shapes, g = 0.)
setPosition!(link1,q = Quaternion(RotX(pi / 2)))
setPosition!(link1,link2,p1 = vert12,p2 = vert21)
# setVelocity!(link1,ω=[2.;0;0])
# setVelocity!(link1,link2,p1=vert12,p2=vert21)

# setForce!([[0.1];nothing],joint1to2,mech)
