using RigidBodyDynamics
using AtlasRobot

# path = "src/util/atlas_simple.urdf"
# path = "src/util/twoTwoBarDiffLength.urdf"
path = "src/util/pendulum.urdf"
doublependulum = RigidBodyDynamics.parse_urdf(path)

state = MechanismState(doublependulum)

ts, qs, vs = simulate(state, 5., Δt = 1e-3);

using MeshCatMechanisms
using Blink

mvis = MechanismVisualizer(doublependulum, URDFVisuals(path));
open(mvis, Blink.Window())
MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 0.001);


