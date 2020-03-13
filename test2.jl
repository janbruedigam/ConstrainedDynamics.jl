using Rotations

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# path = "src/util/atlas_simple.urdf"
# path = "src/util/twoTwoBarDiffLength.urdf"
path = "src/util/pendulum.urdf"
mech = Mechanism(path)
# link1 = mech.bodies[1]
# origin = mech.origin
# joint1 = mech.eqconstraints[2]
# vertices = joint1.constraints[1].vertices


# setPosition!(mech,origin,link1,p1=vertices[1],p2=vertices[2],Δq=Quaternion(RotX(0.2)))

simulate!(mech, save=true)
visualize!(mech)
