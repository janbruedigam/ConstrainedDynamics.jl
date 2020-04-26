using ConstrainedDynamics


path = "examples/examples_files/doublependulum.urdf"
mech = Mechanism(path)

simulate!(mech, save = true)
visualize!(mech)
