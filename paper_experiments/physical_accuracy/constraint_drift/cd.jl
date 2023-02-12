using ConstrainedDynamics
using ConstrainedDynamicsVis
using LinearAlgebra

# Parameters
joint_axis = [1.;0.;0.]
m = 1. # mass
l1 = 1.0 # length1
l2 = sqrt(2)/2 # length2
x,y = .1,.1 # size of link


 # joint connection points
p1 = [0.;0.;l1/2]
p2 = -p1
p3 = [0.;0.;l2/2]
p4 = -p3

# Initial orientation
ϕ1, ϕ2, ϕ3 = 2.842889445268244, -1.209429202890086, pi
q1, q2, q3 = RotX(ϕ1), RotX(ϕ2), RotX(ϕ3)

# Links
origin = Origin()
box1 = Box(x,y,l1,l1, color=RGBA(0,0.3,0.8))
box2 = Box(x,y,l2,l2, color=RGBA(0,0.3,0.8))
box3 = Box(x,y,l1,l1, color=RGBA(0,0.3,0.8))

joint0to13 = EqualityConstraint(Revolute(origin,box1,joint_axis;p2=p1),Revolute(origin,box3,joint_axis;p1=[0;1;0],p2=p1))
joint1to2 = EqualityConstraint(Revolute(box1,box2,joint_axis;p1=p2,p2=p3))
joint2to3 = EqualityConstraint(Revolute(box2,box3,joint_axis;p1=p4,p2=p2))

mech = Mechanism(origin,[box1;box2;box3],[joint0to13;joint1to2;joint2to3])

setPosition!(origin, box1; p2 = p1, Δq = q1)
setPosition!(origin, box3; p1=[0;1;0],p2 = p1, Δq = q3)
setPosition!(box1, box2; p1=p2, p2 = p3, Δq = q2)


drift = zeros(60000)

function controller!(mechanism, k)
    for eqc in mechanism.eqconstraints
        typeof(eqc) <: EqualityConstraint{T,N,Nc,Cs} where {T,N,Nc,Cs<:Tuple{<:Friction}} && continue
        drift[k] += norm(ConstrainedDynamics.gc(mechanism, eqc))
    end
end

storage = simulate!(mech, 600, controller!, record=true)
# visualize(mech, storage)