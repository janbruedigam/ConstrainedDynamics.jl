using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

Δt=0.01
g = -9.81
m = 1.0
l = 1.0
r = 0.01
rod = Cylinder(r, l, m)

p2 = [0.0;0.0;l / 2] # joint connection point

# Links
origin = Origin{Float64}()
link1 = Body(rod)

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, zeros(3), p2, joint_axis))

links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [rod]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=g,Δt=Δt)

T0 = 2*pi*sqrt((m*l^2/3)/(-g*l/2))
t0 = 1

q1 = Quaternion(RotX(π / 2))
setPosition!(mech,origin,link1,p2 = p2,Δq = q1)
simulate!(mech,save = true)
T = T0*1.18
Δ = Int(round(T/Δt))
traj = getindex.(mech.storage.x[1],2)
diff = traj[1:Δ]-traj[Δ+1:2*Δ]
@test maximum(diff)<0.1

q1 = Quaternion(RotX(π / 3))
setPosition!(mech,origin,link1,p2 = p2,Δq = q1)
simulate!(mech,save = true)
T = T0*1.073
Δ = Int(round(T/Δt))
traj = getindex.(mech.storage.x[1],2)
diff = traj[1:Δ]-traj[Δ+1:2*Δ]
@test maximum(diff)<0.1

q1 = Quaternion(RotX(π / 6))
setPosition!(mech,origin,link1,p2 = p2,Δq = q1)
simulate!(mech,save = true)
T = T0*1.017
Δ = Int(round(T/Δt))
traj = getindex.(mech.storage.x[1],2)
diff = traj[1:Δ]-traj[Δ+1:2*Δ]
@test maximum(diff)<0.1



