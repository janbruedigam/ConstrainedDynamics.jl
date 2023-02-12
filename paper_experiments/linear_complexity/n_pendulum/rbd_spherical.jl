using RigidBodyDynamics
using BenchmarkTools
using LinearAlgebra
using StaticArrays


g = -9.81

l = 1.
m = l
I = diagm([0.0845833;0.0845833;0.00125])

# Function to reset joint angles for each run
function reset!(state,jointies,N)
    q =  [0.9238795325112867; 0.3826834323650897; 0.0; 0.0]
    set_configuration!(state, jointies[1], q)
    q =  [1.; 0.0; 0.0; 0.0]
    for i = 2:N
        @eval begin
            set_configuration!($state, $jointies[$i], $q)
        end
    end
    setdirty!(state)
end

timing = zeros(100)

for N = 1:100
    display(N)
    origin = RigidBody{Float64}("origin")
    pendulum = Mechanism(origin; gravity = SVector(0., 0., g))

    joint1 = Joint("joint1", QuaternionSpherical{Float64}())
    inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l/2),moment_about_com=I,mass=m)
    link1 = RigidBody(inertia1)
    before_joint1_to_origin = one(Transform3D,frame_before(joint1), default_frame(origin))
    attach!(pendulum, origin, link1, joint1,joint_pose = before_joint1_to_origin)

    jointies = [joint1]
    links = [link1]

    for i = 2:N
        @eval begin
            $(Symbol("joint",i)) = Joint("joint"*string($i), QuaternionSpherical{Float64}())
            $(Symbol("inertia",i)) = SpatialInertia(frame_after($(Symbol("joint",i))),com=SVector(0, 0, -l/2),moment_about_com=I,mass=m)
            $(Symbol("link",i)) = RigidBody($(Symbol("inertia",i)))
            $(Symbol("before_joint",i,"_to_after_joint",i-1)) = Transform3D(frame_before($(Symbol("joint",i))), frame_after($jointies[$i-1]), SVector(0, 0., -l))
            attach!($pendulum, $links[$i-1], $(Symbol("link",i)), $(Symbol("joint",i)),joint_pose = $(Symbol("before_joint",i,"_to_after_joint",i-1)))

            push!($links,$(Symbol("link",i)))
            push!($jointies, $(Symbol("joint",i)))
        end
    end

    state = MechanismState(pendulum)

    bn = @benchmarkable simulate($state, 10., Δt = 0.01) setup=(reset!($state,$jointies,$N))
    timing[N] = BenchmarkTools.minimum(run(bn,samples=100,seconds=100)).time
end