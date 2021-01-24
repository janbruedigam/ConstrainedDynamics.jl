include("examples/doublependulum_disconnection.jl")

dq = minimalCoordinates(mech)
@test true
setPosition!(mech,dq)
@test true

dF = ConstrainedDynamics.UnitDict(dq.keys,[[0.0] for i=1:2])
setForce!(mech,dF)
@test true