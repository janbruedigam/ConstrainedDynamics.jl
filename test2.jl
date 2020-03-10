!(@isdefined MaximalCoordinateDynamics) && include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics


path = "src/util/twoTwoBarDiffLength.urdf"
mech = Mechanism(path)

simulate!(mech, save=true)
visualize!(mech)