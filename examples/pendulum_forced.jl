using StaticArrays
using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = 0
q1 = Quaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Body(box)

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, zeros(3), p2, joint_axis))

links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes, tend = 10.)
setPosition!(mech,origin,link1,p2 = p2,Δq = q1)

function controller!(mechanism, t)
    τ = SVector{3,Float64}(1, 0, 0) * cos(0.5 * t * 2pi)
    setForce!(mechanism, mechanism.bodies[1], τ = τ)
end


simulate!(mech,controller!,save = true)
visualize!(mech)
