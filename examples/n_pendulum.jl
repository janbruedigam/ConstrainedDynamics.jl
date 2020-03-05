using Rotations
using Plots

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("..", "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
ex = [1.;0.;0.]

h = 1.
r = .05
b1 = Cylinder(r,h,h,color=RGBA(1.,0.,0.))

vert11 = [0.;0.;h/2]
vert12 = -vert11

# Initial orientation
phi = pi/4
q1 = Quaternion(RotX(phi))

# Links
N = 20

origin = Origin{Float64}()
link1 = Body(b1)

links = [link1]

for i=2:N
    @eval begin
        $(Symbol("link",i)) = Body(b1)
        push!(links,$(Symbol("link",i)))
    end
end

# Constraints
jointb1 = EqualityConstraint(Revolute(origin,link1,zeros(3),vert11,ex))

constraints = [jointb1]

for i=2:N
    @eval begin
        $(Symbol("joint",i-1,i)) = EqualityConstraint(Revolute($(Symbol("link",i-1)),$(Symbol("link",i)),vert12,vert11,ex))
        push!(constraints,$(Symbol("joint",i-1,i)))
    end
end

shapes = [b1]

mech = Mechanism(origin,links, constraints;tend=10.,dt=0.01, shapes=shapes)
setPosition!(mech,origin,link1,p2=vert11,Δq=q1)
previd = link1.id
for body in Iterators.drop(mech.bodies,1)
    global previd
    setPosition!(mech,MaximalCoordinateDynamics.getbody(mech,previd),body,p1=vert12,p2=vert11)
    previd = body.id
end

simulate!(mech,save=true)
visualize!(mech)
