using Rotations
using Plots: RGBA

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("..", "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0
l2 = 1.0#sqrt(2)/2
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))
b2 = Box(x,y,l2,l2,color=RGBA(1.,1.,0.))
b3 = Box(x,y,l1,l1,color=RGBA(1.,1.,0))
b4 = Box(x,y,l2,l2,color=RGBA(1.,1.,0))

vert11 = [0.;0.;l1/2]
vert12 = -vert11

vert21 = [0.;0.;l2/2]
vert22 = -vert21

# Initial orientation
offset1 = pi/4
offset2 = pi/2
phi1 = pi/8
q1 = Quaternion(RotX(phi1))
qoff1 = Quaternion(RotX(offset1))
qoff2 = Quaternion(RotX(offset2))

N = 10

# Links
origin = Origin{Float64}()

link1 = Body(b1)
link2 = Body(b2)
link3 = Body(b3)
link4 = Body(b4)


links = [link1;link2;link3;link4]

for i=2:N-1
    @eval begin
        $(Symbol("link",(i-1)*4+1)) = Body(b1)
        push!($links,$(Symbol("link",(i-1)*4+1)))

        $(Symbol("link",(i-1)*4+2)) = Body(b2)
        push!($links,$(Symbol("link",(i-1)*4+2)))

        $(Symbol("link",(i-1)*4+3)) = Body(b3)
        push!($links,$(Symbol("link",(i-1)*4+3)))

        $(Symbol("link",(i-1)*4+4)) = Body(b4)
        push!($links,$(Symbol("link",(i-1)*4+4)))
    end
end

if N>1
    for i=N
        @eval begin
            $(Symbol("link",(i-1)*4+1)) = Body(b1)
            push!($links,$(Symbol("link",(i-1)*4+1)))

            $(Symbol("link",(i-1)*4+2)) = Body(b2)
            push!($links,$(Symbol("link",(i-1)*4+2)))

            $(Symbol("link",(i-1)*4+3)) = Body(b3)
            push!($links,$(Symbol("link",(i-1)*4+3)))

            $(Symbol("link",(i-1)*4+4)) = Body(b4)
            push!($links,$(Symbol("link",(i-1)*4+4)))
        end
    end
end

# Constraints
joint0to1 = EqualityConstraint(Revolute(origin,link1,zeros(3),vert11,ex))
joint1to23 = EqualityConstraint(Revolute(link1,link2,vert12,vert21,ex),Cylindrical(link1,link3,vert11,vert11,ex))
joint3to4 = EqualityConstraint(Revolute(link3,link4,vert12,vert21,ex))
joint2to4 = EqualityConstraint(Revolute(link2,link4,vert22,vert22,ex))

constraints = [joint0to1; joint1to23; joint3to4; joint2to4]

for i=2:N
    @eval begin
        $(Symbol("joint",(i-1)*4,"to",(i-1)*4+1)) = EqualityConstraint(Revolute($links[($i-1)*4],$links[($i-1)*4+1],vert22,vert11,ex))
        push!($constraints,$(Symbol("joint",(i-1)*4,"to",(i-1)*4+1)))

        $(Symbol("joint",(i-1)*4+1,"to",(i-1)*4+2,(i-1)*4+3)) = EqualityConstraint(Revolute($links[($i-1)*4+1],$links[($i-1)*4+2],vert12,vert21,ex),Cylindrical($links[($i-1)*4+1],$links[($i-1)*4+3],vert11,vert11,ex))
        push!($constraints,$(Symbol("joint",(i-1)*4+1,"to",(i-1)*4+2,(i-1)*4+3)))

        $(Symbol("joint",(i-1)*4+3,"to",(i-1)*4+4)) = EqualityConstraint(Revolute($links[($i-1)*4+3],$links[($i-1)*4+4],vert12,vert21,ex))
        push!($constraints,$(Symbol("joint",(i-1)*4+3,"to",(i-1)*4+4)))

        $(Symbol("joint",(i-1)*4+2,"to",(i-1)*4+4)) = EqualityConstraint(Revolute($links[($i-1)*4+2],$links[($i-1)*4+4],vert22,vert22,ex))
        push!($constraints,$(Symbol("joint",(i-1)*4+2,"to",(i-1)*4+4)))
    end
end

shapes = [b1,b2,b3,b4]

mech = Mechanism(origin,links, constraints, shapes=shapes)
if N>1
    setPosition!(mech,origin,link1,p2=vert11,Δq=q1*qoff1)
    setPosition!(mech,link1,link2,p1=vert12,p2=vert21,Δq=inv(q1)*inv(q1))
    setPosition!(mech,link1,link3,p1=vert11,p2=vert11,Δq=inv(q1)*inv(q1))
    setPosition!(mech,link3,link4,p1=vert12,p2=vert21,Δq=q1*q1)
else
    setPosition!(mech,origin,link1,p2=vert11,Δq=q1*qoff2)
    setPosition!(mech,link1,link2,p1=vert12,p2=vert21,Δq=inv(q1)*inv(q1))
    setPosition!(mech,link1,link3,p1=vert11,p2=vert11,Δq=inv(q1)*inv(q1))
    setPosition!(mech,link3,link4,p1=vert12,p2=vert21,Δq=q1*q1)
end

for i=2:N-1
    @eval begin
        setPosition!(mech,$links[($i-1)*4],$(Symbol("link",(i-1)*4+1)),p1=vert22,p2=vert11)
        setPosition!(mech,$(Symbol("link",(i-1)*4+1)),$(Symbol("link",(i-1)*4+2)),p1=vert12,p2=vert21,Δq=inv(q1)*inv(q1))
        setPosition!(mech,$(Symbol("link",(i-1)*4+1)),$(Symbol("link",(i-1)*4+3)),p1=vert11,p2=vert11,Δq=inv(q1)*inv(q1))
        setPosition!(mech,$(Symbol("link",(i-1)*4+3)),$(Symbol("link",(i-1)*4+4)),p1=vert12,p2=vert21,Δq=q1*q1)    
    end
end

if N>1
    for i=N
        @eval begin
            setPosition!(mech,$links[($i-1)*4],$(Symbol("link",(i-1)*4+1)),p1=vert22,p2=vert11,Δq=inv(qoff1)*qoff2)
            setPosition!(mech,$(Symbol("link",(i-1)*4+1)),$(Symbol("link",(i-1)*4+2)),p1=vert12,p2=vert21,Δq=inv(q1)*inv(q1))
            setPosition!(mech,$(Symbol("link",(i-1)*4+1)),$(Symbol("link",(i-1)*4+3)),p1=vert11,p2=vert11,Δq=inv(q1)*inv(q1))
            setPosition!(mech,$(Symbol("link",(i-1)*4+3)),$(Symbol("link",(i-1)*4+4)),p1=vert12,p2=vert21,Δq=q1*q1)       
        end
    end
end


simulate!(mech,save=true)
visualize!(mech)
