!(@isdefined MaximalCoordinateDynamics) && include(joinpath(pwd(), "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

path = "examples/examples_files/doublependulum.urdf"
mech = Mechanism(path)

simulate!(mech, save=true)
visualize!(mech)
