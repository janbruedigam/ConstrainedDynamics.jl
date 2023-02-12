using ConstrainedDynamics
using Statistics

# Parameters
joint_axis = [1.;0.;0.]
m = 1. # mass
l = 1. # length
d = .1 # length of one side of a cube

p1 = [0.;0.;l] # joint connection point


# Initial orientation
ϕ1, ϕ2 = π/2, 0.
q1, q2 = Quaternion(RotX(ϕ1)), Quaternion(RotX(ϕ2))

# Links
origin = Origin{Float64}()

box1 = Box(d,d,d,m,color=RGBA(0,0.3,0.8))
box2 = Box(d,d,d,m,color=RGBA(0,0.3,0.8))


# Constraints
joint0to1 = EqualityConstraint(Revolute(origin,box1,joint_axis;p2=p1))
joint1to2 = EqualityConstraint(Revolute(box1,box2,joint_axis;p2=p1))

# Mechanism
mech = Mechanism(origin, [box1;box2], [joint0to1;joint1to2])

setPosition!(origin, box1; p2 = p1, Δq = q1)
setPosition!(box1, box2; p2 = p1)


energy = zeros(360000)

function controller!(mechanism, k)
    z1 = box1.state.xc[3]
    v1 = box1.state.vc
    ω1 = box1.state.ωc
    T1 = 0.5*v1'*v1*box1.m + 0.5*ω1'*box1.J*ω1
    V1 = box1.m*9.81*z1 

    z2 = box2.state.xc[3]
    v2 = box2.state.vc
    ω2 = box2.state.ωc
    T2 = 0.5*v2'*v2*box2.m + 0.5*ω2'*box2.J*ω2
    V2 = box2.m*9.81*z2
    
    energy[k] = T1+V1+T2+V2
end

storage = simulate!(mech, 3600, controller!, record=true)
visualize(mech,storage)
