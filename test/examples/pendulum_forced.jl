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
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))

links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes)
setPosition!(origin,link1,p2 = p2,Δq = q1)

jointid = constraints[1].id
function pendulum_forced_control!(mechanism, k)
    τ = SVector{1,Float64}(cos(0.5 * k*0.01 * 2pi))
    setForce!(mechanism, geteqconstraint(mechanism,jointid), τ)
    return 
end
