using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedDynamics: g, solve_Eqc,∂g∂posb


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = π / 2
q1 = UnitQuaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Body(box)

# Constraints
#t,r= Revolute(origin, link1, joint_axis; p2=p2)
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))
#println(∂g∂posb(t[1],link1.state.xc,link1.state.qc))
#X,Q=∂g∂posb(r[1],link1.state.xc,link1.state.qc)
#println(X,r[1].V12*Q)
links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes)
setPosition!(origin,link1,p2=[5;6;7],Δq = UnitQuaternion(RotX(0.1)))

solve_Eqc(mech,1e-10,30) 

storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)
