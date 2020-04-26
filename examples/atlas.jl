using ConstrainedDynamics

path = "examples/examples_files/atlas_simple.urdf"
mech = Mechanism(path, floating=true, g = -.5)

simulate!(mech, save = true)
visualize!(mech)
