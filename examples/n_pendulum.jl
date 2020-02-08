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
vert1 = [[vert11];[vert12]]

# Initial orientation
phi = pi/4
q1 = Quaternion(RotX(phi))

# Links
N = 20

origin = Origin{Float64}()

link1 = Body(b1)
setInit!(origin,link1,zeros(3),vert11,q=q1)

links = [link1]

for i=2:N
    @eval begin
        $(Symbol("link",i)) = Body(b1)
        setInit!($(Symbol("link",i-1)),$(Symbol("link",i)),vert12,vert11,q=q1)
        push!(links,$(Symbol("link",i)))
    end
end

# Constraints
jointb1 = Constraint(Revolute(origin,link1,zeros(3),vert11,ex))

constraints = [jointb1]

for i=2:N
    @eval begin
        $(Symbol("joint",i-1,i)) = Constraint(Revolute($(Symbol("link",i-1)),$(Symbol("link",i)),vert12,vert11,ex))
        push!(constraints,$(Symbol("joint",i-1,i)))
    end
end

shapes = [b1]

mech = Mechanism(origin,links, constraints;tend=10.,dt=0.01)

simulate!(mech,save=true)
MaximalCoordinateDynamics.visualize(mech,shapes)
