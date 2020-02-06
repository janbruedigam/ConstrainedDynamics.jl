using Rotations
using BenchmarkTools

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
joint_axis = [1.;0.;0.]
m = 1.
l = 1.0
x,y = .1,.1
box = Box(x,y,l,m)

# joint connection points
p1 = [0.;0.;l/2]
p2 = -p1


# Initial orientation
ϕ1, ϕ2, ϕ3, ϕ4 = π/8,-π/8,-π/8,π/8
q1, q2, q3, q4 = Quaternion(RotX(ϕ1)), Quaternion(RotX(ϕ2)), Quaternion(RotX(ϕ3)), Quaternion(RotX(ϕ4))

offset = π/4
ϕ1, ϕ2, ϕ3, ϕ4 = π/8+offset,-π/8+offset,-π/8+offset,π/8+offset
q5, q6, q7, q8 = Quaternion(RotX(ϕ1)), Quaternion(RotX(ϕ2)), Quaternion(RotX(ϕ3)), Quaternion(RotX(ϕ4))


val = zeros(10)
# Number of segments (-> number of links = 4xN)
for N = 1:10

    # Links
    origin = Origin{Float64}()

    link1 = Link(box)
    link2 = Link(box)
    link3 = Link(box)
    link4 = Link(box)
    if N>1
        setInit!(origin,link1,zeros(3),p1,q=q1)
        setInit!(link1,link2,p2,p1,q=q2)
        setInit!(link1,link3,p1,p1,q=q3)
        setInit!(link3,link4,p2,p1,q=q4)
    else
        setInit!(origin,link1,zeros(3),p1,q=q5)
        setInit!(link1,link2,p2,p1,q=q6)
        setInit!(link1,link3,p1,p1,q=q7)
        setInit!(link3,link4,p2,p1,q=q8)
    end

    links = [link1;link2;link3;link4]

    for i=2:N-1
        @eval begin
            $(Symbol("link",(i-1)*4+1)) = Link(box)
            setInit!($links[($i-1)*4],$(Symbol("link",(i-1)*4+1)),p1,p1,q=q1)
            push!($links,$(Symbol("link",(i-1)*4+1)))

            $(Symbol("link",(i-1)*4+2)) = Link(box)
            setInit!($(Symbol("link",(i-1)*4+1)),$(Symbol("link",(i-1)*4+2)),p2,p1,q=q2)
            push!($links,$(Symbol("link",(i-1)*4+2)))

            $(Symbol("link",(i-1)*4+3)) = Link(box)
            setInit!($(Symbol("link",(i-1)*4+1)),$(Symbol("link",(i-1)*4+3)),p1,p1,q=q3)
            push!($links,$(Symbol("link",(i-1)*4+3)))

            $(Symbol("link",(i-1)*4+4)) = Link(box)
            setInit!($(Symbol("link",(i-1)*4+3)),$(Symbol("link",(i-1)*4+4)),p2,p1,q=q4)
            push!($links,$(Symbol("link",(i-1)*4+4)))
        end
    end

    if N>1
        for i=N
            @eval begin
                $(Symbol("link",(i-1)*4+1)) = Link(box)
                setInit!($links[($i-1)*4],$(Symbol("link",(i-1)*4+1)),p1,p1,q=q5)
                push!($links,$(Symbol("link",(i-1)*4+1)))

                $(Symbol("link",(i-1)*4+2)) = Link(box)
                setInit!($(Symbol("link",(i-1)*4+1)),$(Symbol("link",(i-1)*4+2)),p2,p1,q=q6)
                push!($links,$(Symbol("link",(i-1)*4+2)))

                $(Symbol("link",(i-1)*4+3)) = Link(box)
                setInit!($(Symbol("link",(i-1)*4+1)),$(Symbol("link",(i-1)*4+3)),p1,p1,q=q7)
                push!($links,$(Symbol("link",(i-1)*4+3)))

                $(Symbol("link",(i-1)*4+4)) = Link(box)
                setInit!($(Symbol("link",(i-1)*4+3)),$(Symbol("link",(i-1)*4+4)),p2,p1,q=q8)
                push!($links,$(Symbol("link",(i-1)*4+4)))
            end
        end
    end

    # Constraints
    joint0to1 = Constraint(Socket(origin,link1,zeros(3),p1),Axis(origin,link1,joint_axis))
    joint1to23 = Constraint(Socket(link1,link2,p2,p1),Axis(link1,link2,joint_axis),SocketYZ(link1,link3,p1,p1))
    joint3to4 = Constraint(Socket(link3,link4,p2,p1),Axis(link3,link4,joint_axis))
    joint2to4 = Constraint(Socket(link2,link4,p2,p2),Axis(link2,link4,joint_axis))

    constraints = [joint0to1; joint1to23; joint3to4; joint2to4]

    for i=2:N
        @eval begin
            $(Symbol("joint",(i-1)*4,"to",(i-1)*4+1)) = Constraint(Socket($links[($i-1)*4],$links[($i-1)*4+1],p1,p1),Axis($links[($i-1)*4],$links[($i-1)*4+1],joint_axis))
            push!($constraints,$(Symbol("joint",(i-1)*4,"to",(i-1)*4+1)))

            $(Symbol("joint",(i-1)*4+1,"to",(i-1)*4+2,(i-1)*4+3)) = Constraint(Socket($links[($i-1)*4+1],$links[($i-1)*4+2],p2,p1),Axis($links[($i-1)*4+1],$links[($i-1)*4+2],joint_axis),SocketYZ($links[($i-1)*4+1],$links[($i-1)*4+3],p1,p1))
            push!($constraints,$(Symbol("joint",(i-1)*4+1,"to",(i-1)*4+2,(i-1)*4+3)))

            $(Symbol("joint",(i-1)*4+3,"to",(i-1)*4+4)) = Constraint(Socket($links[($i-1)*4+3],$links[($i-1)*4+4],p2,p1),Axis($links[($i-1)*4+3],$links[($i-1)*4+4],joint_axis))
            push!($constraints,$(Symbol("joint",(i-1)*4+3,"to",(i-1)*4+4)))

            $(Symbol("joint",(i-1)*4+2,"to",(i-1)*4+4)) = Constraint(Socket($links[($i-1)*4+2],$links[($i-1)*4+4],p2,p2),Axis($links[($i-1)*4+2],$links[($i-1)*4+4],joint_axis))
            push!($constraints,$(Symbol("joint",(i-1)*4+2,"to",(i-1)*4+4)))
        end
    end

    qin = [links[i].q[1] for i=1:N*4]
    xin = [links[i].x[1] for i=1:N*4]

    # mechanism
    bot = Robot(origin,links, constraints)

    t = @benchmarkable simulate_reset!($bot,$xin,$qin)
    val[N] = BenchmarkTools.minimum(run(t,samples=100,seconds=100)).time
end
