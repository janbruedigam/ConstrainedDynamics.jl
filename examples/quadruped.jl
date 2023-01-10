using ConstrainedDynamics
using ConstrainedDynamicsVis
using StaticArrays
using ConstrainedDynamics: g
using Plots


path = "examples/examples_files/quadruped_simple.urdf"
mech = Mechanism(path, floating=false, g = 0)
origin = mech.origin
bodies = mech.bodies.values
eqcs = mech.eqconstraints.values

fricsandineqsFR = Friction(getbody(mech,"FR_calf"), [0;0;1.0], 0.2; p = [0.0;0;-0.1])
fricsandineqsFL = Friction(getbody(mech,"FL_calf"), [0;0;1.0], 0.2; p = [0.0;0;-0.1])
fricsandineqsRR = Friction(getbody(mech,"RR_calf"), [0;0;1.0], 0.2; p = [0.0;0;-0.1])
fricsandineqsRL = Friction(getbody(mech,"RL_calf"), [0;0;1.0], 0.2; p = [0.0;0;-0.1])
fricsFR = fricsandineqsFR[1]
fricsFL = fricsandineqsFL[1]
fricsRR = fricsandineqsRR[1]
fricsRL = fricsandineqsRL[1]
ineqcFR = fricsandineqsFR[2]
ineqcFL = fricsandineqsFL[2]
ineqcRR = fricsandineqsRR[2]
ineqcRL = fricsandineqsRL[2]

frics = [fricsFR;fricsFL;fricsRR;fricsRL]
ineqcs = [ineqcFR;ineqcFL;ineqcRR;ineqcRL]


mech2 = Mechanism(origin, bodies, eqcs, ineqcs, frics, Δt = 0.001)

# setPosition!(mech2, geteqconstraint(mech2,"floating_base"),[0;0;0.232;0.;0.;0.])
setPosition!(mech2, geteqconstraint(mech2,"floating_base"),[-0.007044513654001689, 0.0028495585299927887, 0.21575955182909478,0.00102226, -0.0247454, 0.00297319])
# setPosition!(getbody(mech2,"trunk"),x=[-0.007044513654001689, 0.0028495585299927887, 0.21175955182909478],q=Quaternion(0.999886, 0.00102226, -0.0147454, 0.00297319))


initangle = 0.95

# setPosition!(mech2, geteqconstraint(mech2,"FR_thigh_joint"),[initangle])
# setPosition!(mech2, geteqconstraint(mech2,"FR_calf_joint"),[-2*initangle])

# setPosition!(mech2, geteqconstraint(mech2,"FL_thigh_joint"),[initangle*0.9])
# setPosition!(mech2, geteqconstraint(mech2,"FL_calf_joint"),[-2*initangle])

# setPosition!(mech2, geteqconstraint(mech2,"RR_thigh_joint"),[initangle*0.9])
# setPosition!(mech2, geteqconstraint(mech2,"RR_calf_joint"),[-2*initangle])

# setPosition!(mech2, geteqconstraint(mech2,"RL_thigh_joint"),[initangle])
# setPosition!(mech2, geteqconstraint(mech2,"RL_calf_joint"),[-2*initangle])

setPosition!(mech2, geteqconstraint(mech2,"FR_thigh_joint"),[0.9378205336153784])
setPosition!(mech2, geteqconstraint(mech2,"FR_calf_joint"),[-1.9767045221686401])

setPosition!(mech2, geteqconstraint(mech2,"FL_thigh_joint"),[0.8411651545116884])
setPosition!(mech2, geteqconstraint(mech2,"FL_calf_joint"),[-1.9583154751653324])

setPosition!(mech2, geteqconstraint(mech2,"RR_thigh_joint"),[0.8265886839172527])
setPosition!(mech2, geteqconstraint(mech2,"RR_calf_joint"),[-2.0094918578851146])

setPosition!(mech2, geteqconstraint(mech2,"RL_thigh_joint"),[0.9294618666565416])
setPosition!(mech2, geteqconstraint(mech2,"RL_calf_joint"),[-2.026072700146945])

N = 10000
steps = Base.OneTo(N)
storagequad = Storage(steps,length(bodies))

traj21 = [initangle*(cos(i*0.001*2*pi)*0.1+0.9) for i=1:N]
# traj21 = [[initangle for i=1:1];traj21]
traj31 = [-2*initangle-sin(i*0.001*2*pi)*0.1 for i=1:N]
# traj31 = [[-2*initangle for i=1:10];traj31]

traj22 = [initangle*(cos(i*0.001*2*pi+pi)*0.1+0.9) for i=1:N]
# traj22 = [[initangle*0.9 for i=1:10];traj22]
traj32 = [-2*initangle-sin(i*0.001*2*pi+pi)*0.1 for i=1:N]
# traj32 = [[-2*initangle for i=1:10];traj32]



function singleleg(mechanism, leg, angles)
    j1 = geteqconstraint(mechanism, leg*"_hip_joint")
    j2 = geteqconstraint(mechanism, leg*"_thigh_joint")
    j3 = geteqconstraint(mechanism, leg*"_calf_joint")

    θ1 = minimalCoordinates(mechanism, j1)[1]
    θ2 = minimalCoordinates(mechanism, j2)[1]
    θ3 = minimalCoordinates(mechanism, j3)[1]
    dθ1 = minimalVelocities(mechanism, j1)[1]
    dθ2 = minimalVelocities(mechanism, j2)[1]
    dθ3 = minimalVelocities(mechanism, j3)[1]

    u1 = 100.0*(angles[1]-θ1) + 5.0*(0-dθ1)
    u2 = 80.0*(angles[2]-θ2) + 4.0*(0-dθ2)
    u3 = 60.0*(angles[3]-θ3) + 3.0*(0-dθ3)

    setForce!(mechanism, j1, SA[u1])
    setForce!(mechanism, j2, SA[u2])
    setForce!(mechanism, j3, SA[u3])
end

function controller!(mechanism, k)
    singleleg(mechanism, "FR", SA[0.0;traj21[k];traj31[k]])
    singleleg(mechanism, "FL", SA[0.0;traj22[k];traj32[k]])
    singleleg(mechanism, "RR", SA[0.0;traj22[k];traj32[k]])
    singleleg(mechanism, "RL", SA[0.0;traj21[k];traj31[k]])
end


storagequad = simulate!(mech2, storagequad, controller!, ε = 1e-6)
visualize(mech2, storagequad, showframes = false)

# plot([g(ineqcRR[1].constraints[1],storage.x[2][i],storage.q[2][i])[1] for i=1:N])
# plot([g(ineqcRL[1].constraints[1],storage.x[6][i],storage.q[6][i])[1] for i=1:N])
# plot([g(ineqcFR[1].constraints[1],storage.x[10][i],storage.q[10][i])[1] for i=1:N])
# plot([g(ineqcRR[1].constraints[1],storage.x[13][i],storage.q[13][i])[1] for i=1:N])