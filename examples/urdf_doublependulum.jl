using ConstrainedDynamics


path = "examples/examples_files/doublependulum.urdf"
mech, shapes = Mechanism(path)

storage = simulate!(mech, 10., record = true)
visualize!(mech, storage, shapes)
