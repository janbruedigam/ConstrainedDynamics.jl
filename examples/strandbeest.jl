### Does not work yet!

using ConstrainedDynamics
using ConstrainedDynamicsVis


path = "examples/examples_files/strandbeest.urdf"
mech, shapes = Mechanism(path, floating=true, g = 0)

v = rand.(length.(minimalCoordinates(mech).values))
# d = ConstrainedDynamics.UnitDict(92:92+90,v)
d = ConstrainedDynamics.UnitDict(93:93+91,v)
setPosition!(mech,d)

storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)


