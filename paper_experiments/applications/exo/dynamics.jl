using ConstrainedDynamics
using ConstrainedDynamicsVis
using StaticArrays


path = "experiments/applications/exo/deps/model.urdf"
mech = Mechanism(path, floating=false)
origin = mech.origin
bodies = mech.bodies.values
eqcs = mech.eqconstraints.values

upper_connection = EqualityConstraint(ConstrainedDynamics.PlanarAxis(mech.bodies["sIE_Unit"], mech.bodies["upper_arm"],[0;1;0];p1=[0.0006;-0.1079;-0.1237], spring=100))
lower_connection = EqualityConstraint(ConstrainedDynamics.PlanarAxis(mech.bodies["eFE_Unit"], mech.bodies["lower_arm"],[0;0;1];p1=[0;-0.0693;0.0803],p2=[0;0;0], spring=100))

mech = Mechanism(origin, bodies, [eqcs;upper_connection;lower_connection], g=-9.81, Δt = 0.001)

q = RotX(pi/2) * RotZ(-pi/4)
axisangle = ConstrainedDynamics.rotation_vector(q)
setPosition!(mech, mech.eqconstraints["shoulder"], axisangle, iter=false)
setPosition!(mech, mech.eqconstraints["elbow"], SA[-pi/4], iter=false)
setPosition!(mech, mech.eqconstraints["sFE"], SA[pi/4], iter=false)
setPosition!(mech, mech.eqconstraints["sIE"], SA_F64[], iter=false)
setPosition!(mech, mech.eqconstraints["eFE"], SA[-pi/4], iter=false)

u1 = zeros(10000)
u2 = zeros(10000)

Kh = [100;50]
Kl = Kh/4.75

function controller!(mechanism, k)
    θelbow = minimalCoordinates(mechanism, mech.eqconstraints["elbow"])[1]
    ωelbow = minimalVelocities(mechanism, mech.eqconstraints["elbow"])[1]
    ωeFE = minimalVelocities(mechanism, mech.eqconstraints["eFE"])[1]
    θsFE = minimalCoordinates(mechanism, mech.eqconstraints["sFE"])[1]
    θeFE = minimalCoordinates(mechanism, mech.eqconstraints["eFE"])[1]
    ωsFE = minimalVelocities(mechanism, mech.eqconstraints["sFE"])[1]
    ωeFE = minimalVelocities(mechanism, mech.eqconstraints["eFE"])[1]

    uelbow = 10*(-pi/4-θelbow) + 2*(0-ωelbow) # spastic
    # uelbow = 0 # healthy

    θsFE_des = pi/4 + pi/4*sin(0.001*k)
    θeFE_des = -pi/4 + pi/4*sin(0.001*k)
    u1[k] = usFE = Kl[1]*(θsFE_des-θsFE) + Kl[1]/10*(0-ωsFE)
    u2[k] = ueFE = Kl[2]*(θeFE_des-θeFE) + Kl[2]/10*(0-ωeFE)
    # u1[k] = usFE = Kh[1]*(θsFE_des-θsFE) + Kh[1]/10*(0-ωsFE)
    # u2[k] = ueFE = Kh[2]*(θeFE_des-θeFE) + Kh[2]/10*(0-ωeFE)

    setForce!(mechanism, mech.eqconstraints["elbow"], [uelbow])
    setForce!(mechanism, mech.eqconstraints["sFE"], [usFE])
    setForce!(mechanism, mech.eqconstraints["eFE"], [ueFE])
end

storage = simulate!(mech, 10, controller!, record = true)
visualize(mech, storage)