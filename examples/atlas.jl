using ConstrainedDynamics


path = "examples/examples_files/atlas_simple.urdf"
mech, shapes = Mechanism(path, floating=true, g = -.5)

storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)
