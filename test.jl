using Rotations

(@isdefined MaximalCoordinateDynamics) ? nothing : include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
l1 = 1.0
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))

x1 = [1.;2.;3.]
x2 = [-1.;2.;-3.]
q1 = Quaternion(RotX(1.))
q2 = Quaternion(RotX(2.))

pa = [1.;2.;3.]
pb = [4.;5.;6.]

# Bodies
origin = Origin{Float64}()

link1 = Body(b1)
setInit!(origin,link1,x1,zeros(3),q=q1)

link2 = Body(b1)
setInit!(origin,link2,x2,zeros(3),q=q2)


# Constraints

shapes = [b1]

oc1 = Constraint(OriginConnection(origin,link1))
oc2 = Constraint(OriginConnection(origin,link2))

joint1 = Constraint(Socket(link1,link2,pa,pb))

bot = Mechanism(origin,[link1;link2],[oc1;oc2;joint1])

# simulate!(bot,save=true,debug=false)
# MaximalCoordinateDynamics.visualize(bot,shapes)
