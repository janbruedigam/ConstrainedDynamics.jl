!(@isdefined MaximalCoordinateDynamics) && include(joinpath(pwd(), "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

path = "examples/examples_files/atlas_simple.urdf"
mech = Mechanism(path, g = -.5)

simulate!(mech, save = true)
visualize!(mech)
