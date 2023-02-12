using ConstrainedDynamics
using ConstrainedDynamicsVis
using StaticArrays
using Random
using LinearAlgebra


path = "experiments/applications/walker/deps/quadruped_simple.urdf"
N = 10000
cf = 0.2
contact_normal = [0;0;1]
contact_point = [0;0;-0.1]

mech = Mechanism(path, floating=false)
origin = mech.origin
bodies = mech.bodies.values
eqcs = mech.eqconstraints.values

fricsandineqsFR = Friction(getbody(mech,"FR_calf"), contact_normal, cf; p = contact_point)
fricsandineqsFL = Friction(getbody(mech,"FL_calf"), contact_normal, cf; p = contact_point)
fricsandineqsRR = Friction(getbody(mech,"RR_calf"), contact_normal, cf; p = contact_point)
fricsandineqsRL = Friction(getbody(mech,"RL_calf"), contact_normal, cf; p = contact_point)

frics = [fricsandineqsFR[1];fricsandineqsFL[1];fricsandineqsRR[1];fricsandineqsRL[1]]
ineqcs = [fricsandineqsFR[2];fricsandineqsFL[2];fricsandineqsRR[2];fricsandineqsRL[2]]

mech = Mechanism(origin, bodies, eqcs, ineqcs, frics, Δt = 0.001)

storagequad = Storage(N,length(bodies))

calf_body = getbody(mech, "FR_calf")
trunk_body = getbody(mech, "trunk")
floating_constraint = geteqconstraint(mech, "floating_base")

legjoints = [
    geteqconstraint(mech, "FR_hip_joint")
    geteqconstraint(mech, "FR_thigh_joint")
    geteqconstraint(mech, "FR_calf_joint")
    geteqconstraint(mech, "FL_hip_joint")
    geteqconstraint(mech, "FL_thigh_joint")
    geteqconstraint(mech, "FL_calf_joint")
    geteqconstraint(mech, "RR_hip_joint")
    geteqconstraint(mech, "RR_thigh_joint")
    geteqconstraint(mech, "RR_calf_joint")
    geteqconstraint(mech, "RL_hip_joint")
    geteqconstraint(mech, "RL_thigh_joint")
    geteqconstraint(mech, "RL_calf_joint")
]
