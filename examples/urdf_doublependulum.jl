!(@isdefined MaximalCoordinateDynamics) && include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

path = "examples/examples_files/doublependulum.urdf"
mech = Mechanism(path)

simulate!(mech, save=true)
visualize!(mech)
