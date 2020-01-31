using LinearAlgebra
using RigidBodyDynamics
using StaticArrays
using BenchmarkTools
using Rotations
using LinearAlgebra
using Blink

using MeshCatMechanisms

g = -9.81

l1 = 1.
m1 = l1
I1 = diagm([0.0845833;0.0845833;0.00125])

# axis = SVector(1., 0., 0.)

function test(state,jointies,N)
    ang = 0.2
    qu = Quat(RotX(ang))
    q1 = [qu.w;qu.x;qu.y;qu.z]
    qu = Quat(RotX(pi/2-ang))
    q2 = [qu.w;qu.x;qu.y;qu.z]
    qu = Quat(RotX(pi/2+ang))
    q3 = [qu.w;qu.x;qu.y;qu.z]

    set_configuration!(state, jointies[1], q1)
    for i = 2:floor(Int64,N/2)
        set_configuration!(state, jointies[i], [1.;0;0;0])
    end
    set_configuration!(state, jointies[ceil(Int64,N/2)], q2)
    set_configuration!(state, jointies[ceil(Int64,N/2)+1], q3)
    for i = ceil(Int64,N/2)+2:N-1
        set_configuration!(state, jointies[i], [1.;0;0;0])
    end
    set_configuration!(state, jointies[N], q1)

    setdirty!(state)
    simulate(state, 10., Δt = 0.01)
end

# val2 = zeros(101)

for N = [3;5;7;9;11;21]
    world = RigidBody{Float64}("world")
    npend = Mechanism(world; gravity = SVector(0., 0., g))

    joint1 = Joint("joint1", QuaternionSpherical{Float64}())
    inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l1/2),moment_about_com=I1,mass=m1)
    link1 = RigidBody(inertia1)
    before_joint1_to_world = one(Transform3D,frame_before(joint1), default_frame(world))
    attach!(npend, world, link1, joint1,joint_pose = before_joint1_to_world)

    jointies = [joint1]
    links = [link1]

    for i = 2:N-1
        @eval begin
            $(Symbol("joint",i)) = Joint("joint"*string($i), QuaternionSpherical{Float64}())
            $(Symbol("inertia",i)) = SpatialInertia(frame_after($(Symbol("joint",i))),com=SVector(0, 0, -l1/2),moment_about_com=I1,mass=m1)
            $(Symbol("link",i)) = RigidBody($(Symbol("inertia",i)))
            $(Symbol("before_joint",i,"_to_after_joint",i-1)) = Transform3D(frame_before($(Symbol("joint",i))), frame_after($jointies[$i-1]), SVector(0, 0., -l1))
            attach!($npend, $links[$i-1], $(Symbol("link",i)), $(Symbol("joint",i)),joint_pose = $(Symbol("before_joint",i,"_to_after_joint",i-1)))

            push!($links,$(Symbol("link",i)))
            push!($jointies, $(Symbol("joint",i)))
        end
    end

    @eval begin
        $(Symbol("joint",N)) = Joint("joint"*string($N), QuaternionSpherical{Float64}())
        $(Symbol("inertia",N)) = SpatialInertia(frame_after($(Symbol("joint",N))),com=SVector(0, 0, -l1/2),moment_about_com=I1,mass=m1)
        $(Symbol("link",N)) = RigidBody($(Symbol("inertia",N)))
        $(Symbol("before_joint",N,"_to_world")) = Transform3D(frame_before($(Symbol("joint",N))), default_frame($world), SVector(0., 1.,0.))
        attach!($npend, $world, $(Symbol("link",N)), $(Symbol("joint",N)),joint_pose = $(Symbol("before_joint",N,"_to_world")))
        push!($links,$(Symbol("link",N)))
        push!($jointies,$(Symbol("joint",N)))
    end

    @eval begin
        $(Symbol("joint",N+1)) = Joint("joint"*string($N+1), QuaternionSpherical{Float64}())
        $(Symbol("before_joint",N+1,"_to_joint",N-1)) = Transform3D(frame_before($(Symbol("joint",N+1))), frame_after($jointies[$N-1]), SVector(0, 0., -l1))
        $(Symbol("joint",N,"_to_joint",N+1)) = Transform3D(frame_after($jointies[$N]), frame_after($(Symbol("joint",N+1))), SVector(0, 0., l1))
        attach!($npend, $links[$N-1], $links[$N], $(Symbol("joint",N+1)),joint_pose = $(Symbol("before_joint",N+1,"_to_joint",N-1)), successor_pose = $(Symbol("joint",N,"_to_joint",N+1)))
        push!($jointies,$(Symbol("joint",N+1)))
    end

    state = MechanismState(npend)

    t = @benchmarkable test($state,$jointies,$N)
    val2[N] = BenchmarkTools.minimum(run(t,samples=100,seconds=100)).time
end
