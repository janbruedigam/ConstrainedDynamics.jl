using ConstrainedDynamics
using LinearAlgebra
using StaticArrays


# Parameters
h = 0.06
r = 0.08
b1 = Cylinder(r*0.4, h, h*(r*0.4)^2*500, color = RGBA(0.8, 0.2, 0.))
b2 = Cylinder(r*0.8, h, h*(r*0.8)^2*500, color = RGBA(0.8, 0.2, 0.))
b3 = Cylinder(r, h, h*r^2*500, color = RGBA(0.8, 0.2, 0.))

# Links
origin = Origin{Float64}()
link1 = Body(b1)
link2 = Body(b2)
link3 = Body(b3)
link4 = Body(b2)
link5 = Body(b1)


# Constraints
joint1 = EqualityConstraint(Fixed(link2, link1; p1=-[0;0;0.03], p2=[0;0;0.03]))
joint2 = EqualityConstraint(Fixed(link3, link2; p1=-[0;0;0.03], p2=[0;0;0.03]))
joint3 = EqualityConstraint(OriginConnection(origin, link3))
joint4 = EqualityConstraint(Fixed(link3, link4; p1=[0;0;0.03], p2=-[0;0;0.03]))
joint5 = EqualityConstraint(Fixed(link4, link5; p1=[0;0;0.03], p2=-[0;0;0.03]))
fr1 = InequalityConstraint(Friction(link1, [0;0;1.0], 0.6))
fr2 = InequalityConstraint(Friction(link2, [0;0;1.0], 0.6))
fr3 = InequalityConstraint(Friction(link3, [0;0;1.0], 0.6))
fr4 = InequalityConstraint(Friction(link4, [0;0;1.0], 0.6))
fr5 = InequalityConstraint(Friction(link5, [0;0;1.0], 0.6))

links = [link1;link2;link3;link4;link5]
constraints = [joint1;joint2;joint3;joint4;joint5]
fr = [fr1;fr2;fr3;fr4;fr5]
shapes = [b1;b2;b3]


mech = Mechanism(origin, links, constraints, fr, shapes = shapes)
setPosition!(link3,x = [0.;-10.0;1.5],q=Quaternion(RotX(pi/2)))
setPosition!(link3,link2,Δx = -[0.;0.;0.03])
setPosition!(link2,link1,Δx = -[0.;0.;0.03])
setPosition!(link3,link4,Δx = [0.;0.;0.03])
setPosition!(link4,link5,Δx = [0.;0.;0.03])

spin = 0.35

function fottball_control!(mechanism, k)
    if k<25
        setForce!(link3, F = SA[0.;25.;25.], τ=spin*SA[0.2;0.2;1.])
    elseif k==40
        setForce!(link3, τ = SA[0.;0.3;0.])
    else
        setForce!(link3)
    end
    return
end
