using ConstrainedDynamics

path = "examples/examples_files/atlas_simple.urdf"
mech, shapes = Mechanism(path, floating=false, g = -.5)
