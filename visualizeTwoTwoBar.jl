using RigidBodyDynamics, MeshCatMechanisms, Blink, LinearAlgebra, StaticArrays

#srcdir = dirname(pathof(RigidBodyDynamics))
#urdf = joinpath(srcdir, "..", "test", "urdf", "Acrobot.urdf")
urdf = "twoTwoBarDiffLength.urdf";

grav = -9.81 # gravitational acceleration in z-direction
mechanism = RigidBodyDynamics.parse_urdf(urdf; gravity = SVector(0, 0, grav))

# shoulder, elbow = joints(mechanism)
# function simple_control!(torques::AbstractVector, t, state::MechanismState)
#     torques[velocity_range(state, shoulder)] .= 0 * -1 .* velocity(state, shoulder)
#     torques[velocity_range(state, elbow)] .= 0 * 10 * sin(t)
# end;


jointies = joints(mechanism)
#1 3 4 2


state = MechanismState(mechanism)
#zero_velocity!(state)
set_configuration!(state, jointies[1], pi/2)
set_configuration!(state, jointies[3], -pi/4)
set_configuration!(state, jointies[2], 0)
set_configuration!(state, jointies[4], 3*pi/3)

final_time = 10.
ts, qs, vs = simulate(state, final_time);
ts2, qs2, vs2 = simulate(state, final_time);

for i=1:100000
    qs2[i][1] = trajS[1,Int(ceil(i/100))];
    qs2[i][3] = trajS[2,Int(ceil(i/100))]-trajS[1,Int(ceil(i/100))];
    qs2[i][2] = trajS[3,Int(ceil(i/100))];
    qs2[i][4] = trajS[4,Int(ceil(i/100))]-trajS[3,Int(ceil(i/100))];
end
qs2[100001][1] = qs2[100000][1]
qs2[100001][3] = qs2[100000][3]
qs2[100001][2] = qs2[100000][2]
qs2[100001][4] = qs2[100000][4]

qs2[100002][1] = qs2[100000][1]
qs2[100002][3] = qs2[100000][3]
qs2[100002][2] = qs2[100000][2]
qs2[100002][4] = qs2[100000][4]



# mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf));
mvis2 = MechanismVisualizer(mechanism, URDFVisuals(urdf));
# open(mvis, Blink.Window())
open(mvis2, Blink.Window())
# MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 1.)
MeshCatMechanisms.animate(mvis2, ts2, qs2; realtimerate = 1.)

qsTraj = zeros(4,100002);
for i=1:100002
    qsTraj[1,i] = qs[i][1];
    qsTraj[4,i] = qs[i][2];
    qsTraj[2,i] = qs[i][3];
    qsTraj[3,i] = qs[i][4];
end
