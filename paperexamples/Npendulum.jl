using Rotations
using Plots
using BenchmarkTools

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
h = 1.
r = .05
b1 = Cylinder(r,h,h,color=RGBA(1.,0.,0.))

vert11 = [0.;0.;h/2]
vert12 = -vert11


ex = [1.;0.;0.]



shapes = [b1]

# val = zeros(100)


# for N=[[i for i=1:9];[i for i=10:10:70]]
for N=90:10:100
    ang = pi/4#pi/4 - (rand()*pi/2)
    q1 = Quaternion(RotX(ang))#Quaternion(rand(RotMatrix{3}))
    for i=2:N
        @eval begin
            $(Symbol("q",i)) = Quaternion(RotX($ang))#Quaternion(rand(RotMatrix{3}))
        end
    end

    # Links
    origin = Origin{Float64}()
    link1 = Link(b1)

    xin =[link1.x[1]]
    qin = [q1 for i=1:N]

    links = [link1]

    setInit!(origin,link1,zeros(3),vert11,q=q1)
    for i=2:N
        @eval begin
            $(Symbol("link",i)) = Link(b1)
            push!($links,$(Symbol("link",i)))
            setInit!($links[$i-1],$links[$i],vert12,vert11,q=$(Symbol("q",i)))
            push!($xin,$links[$i].x[1])
        end
    end

    # Constraints
    jointb1 = Constraint(Socket(origin,link1,zeros(3),vert11),Axis(origin,link1,ex))

    constraints = [jointb1]

    for i=2:N
        @eval begin
            $(Symbol("joint",i-1,i)) = Constraint(Socket($links[$i-1],$links[$i],vert12,vert11),Axis($links[$i-1],$links[$i],ex))
            push!($constraints,$(Symbol("joint",i-1,i)))
        end
    end

    bot = Robot(origin,links, constraints;tend=10.,dt=0.01)

    t = @benchmarkable simulate!($bot,$xin,$qin)
    val[N] = BenchmarkTools.median(run(t,samples=100,seconds=100)).time
end
