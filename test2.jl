!(@isdefined MaximalCoordinateDynamics) && include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics


path = "src/util/twoTwoBarDiffLength.urdf"
Mechanism(path)