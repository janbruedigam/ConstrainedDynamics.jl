using Rotations
using Plots

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))

vert11 = [0.;0.;l1/2]
vert12 = -vert11

# Initial orientation
N = 50
phi = 0.5
offset = 0.2
q0 = Quaternion(RotX(0.0+offset))
q1 = Quaternion(RotX(phi+offset))
q2 = Quaternion(RotX(-phi+offset))

# Links
origin = Origin{Float64}()

link1 = Link(b1)
setInit!(origin,link1,zeros(3),vert11,q=q1)
links = [link1]

xvals = [link1.x[1]]
qvals = [[q1];[q0 for i=1:Int64(N/2)-1];[q2];[q1];[q0 for i=1:Int64(N/2)-1];[q2]]

for i=2:Int64(N/2)-1
    @eval begin
        $(Symbol("link",i)) = Link(b1)
        setInit!($(Symbol("link",i-1)),$(Symbol("link",i)),vert12,vert11,q=q0)
        push!(links,$(Symbol("link",i)))
        push!(xvals,$(Symbol("link",i)).x[1])
    end
end

@eval begin
    $(Symbol("link",Int64(N/2))) = Link(b1)
    setInit!($(Symbol("link",Int64(N/2)-1)),$(Symbol("link",Int64(N/2))),vert12,vert11,q=q2)
    push!(links,$(Symbol("link",Int64(N/2))))
    push!(xvals,$(Symbol("link",Int64(N/2))).x[1])
end

@eval begin
    $(Symbol("link",Int64(N/2)+1)) = Link(b1)
    setInit!(origin,$(Symbol("link",Int64(N/2)+1)),zeros(3),vert11,q=q2)
    push!(links,$(Symbol("link",Int64(N/2)+1)))
    push!(xvals,$(Symbol("link",Int64(N/2)+1)).x[1])
end

for i=Int64(N/2)+2:N-1
    @eval begin
        $(Symbol("link",i)) = Link(b1)
        setInit!($(Symbol("link",i-1)),$(Symbol("link",i)),vert12,vert11,q=q0)
        push!(links,$(Symbol("link",i)))
        push!(xvals,$(Symbol("link",i)).x[1])
    end
end

@eval begin
    $(Symbol("link",N)) = Link(b1)
    setInit!($(Symbol("link",N-1)),$(Symbol("link",N)),vert12,vert11,q=q1)
    push!(links,$(Symbol("link",N)))
    push!(xvals,$(Symbol("link",N)).x[1])
end

# Constraints
@eval begin
    $(Symbol("joint0to1",Int64(N/2)+1)) = Constraint(Socket(origin,link1,zeros(3),vert11),Axis(origin,link1,ex),SocketYZ(origin,links[Int64(N/2)+1],zeros(3),vert11))
    joint12 = Constraint(Socket(links[1],links[2],vert12,vert11),Axis(links[1],links[2],ex))
    constraints = [$(Symbol("joint0to1",Int64(N/2)+1));joint12]
end
for i=3:Int64(N/2)
    @eval begin
        $(Symbol("joint",i-1,i)) = Constraint(Socket(links[$i-1],links[$i],vert12,vert11),Axis(links[$i-1],links[$i],ex))
        push!(constraints,$(Symbol("joint",i-1,i)))
    end
end
for i=Int64(N/2)+2:N
    @eval begin
        $(Symbol("joint",i-1,i)) = Constraint(Socket(links[$i-1],links[$i],vert12,vert11),Axis(links[$i-1],links[$i],ex))
        push!(constraints,$(Symbol("joint",i-1,i)))
    end
end
@eval begin
    $(Symbol("joint",Int64(N/2),N)) = Constraint(Socket(links[Int64(N/2)],links[N],vert12,vert12),Axis(links[Int64(N/2)],links[N],ex))
    push!(constraints,$(Symbol("joint",Int64(N/2),N)))
end

shapes = [b1]

bot = Robot(origin,links, constraints,dt=0.01)

drift2 = simulate!(bot,save=true,debug=false)
# FullCordDynamics.visualize(bot,shapes)
