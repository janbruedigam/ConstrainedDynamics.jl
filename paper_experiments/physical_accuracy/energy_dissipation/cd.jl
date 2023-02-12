using ConstrainedDynamics
using ConstrainedDynamicsVis


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = π / 2
q1 = Quaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1)

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2, damper=0.5))

links = [link1]
constraints = [joint_between_origin_and_link1]


mech = Mechanism(origin, links, constraints, Δt = 0.05)
setPosition!(origin,link1,p2 = p2,Δq = q1)

angle = zeros(200)
energy = zeros(200)
function controller!(mechanism, k)
    angle[k] = minimalCoordinates(mech, joint_between_origin_and_link1)[1]

    z1 = link1.state.xc[3]
    v1 = link1.state.vc
    ω1 = link1.state.ωc
    T1 = 0.5*v1'*v1*link1.m + 0.5*ω1'*link1.J*ω1
    V1 = link1.m*9.81*z1 

    energy[k] = T1+V1
end

storage = simulate!(mech, 10., controller!, record = true)
# visualize(mech, storage)