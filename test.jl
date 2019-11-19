include("quaternion.jl")
include("link.jl")
include("constraint.jl")
include("robot.jl")
include("newton.jl")
include("simulation.jl")
include("shapes.jl")

# Parameters
ax = [1.;0.;0.]

l1 = 1.
m1, J1 = box(.1,.1,l1,l1)
vert11 = [0.;0.;l1/2]
vert12 = -vert11
vert1 = [[vert11];[vert12]]

l2 = sqrt(2)/2
m2, J2 = box(.1,.1,l2,l2)
vert21 = [0.;0.;l2/2]
vert22 = -vert21
vert2 = [[vert21];[vert22]]

# Initial orientation
phi1, phi2, phi3, phi4 = pi/2, -pi/4, 0., 3*pi/4.
q1, q2, q3, q4 = Quaternion(RotX(phi1)), Quaternion(RotX(phi2)), Quaternion(RotX(phi3)), Quaternion(RotX(phi4))

# Links
baseLink = Link(Float64)

link1 = Link(m1,J1,vert1)
x1 =  initialPosition(baseLink,link1,q1,[1;1],1)
setInit!(link1,x=x1,q=q1)

link2 = Link(m2,J2,vert2)
x2 = initialPosition(link1,link2,q2,[2;1],1)
setInit!(link2,x=x2,q=q2)

link3 = Link(m1,J1,vert1)
x3 = initialPosition(baseLink,link3,q3,[1;1],1)
setInit!(link3,x=x3,q=q3)

link4 = Link(m2,J2,vert2)
x4 = initialPosition(link3,link4,q4,[2;1],1)
setInit!(link4,x=x4,q=q4)

# Constraints
jointb = [SocketConstraint(baseLink,1); RotationConstraint(baseLink)]
jointb1 = [SocketConstraint(baseLink,link1,[1;1]); RotationConstraint(baseLink,link1,ax)]
joint12 = [SocketConstraint(link1,link2,[2;1]); RotationConstraint(link1,link2,ax)]
jointb3 = [SocketConstraint(baseLink,link3,[1;1]); RotationConstraint(baseLink,link3,ax)]
joint34 = [SocketConstraint(link3,link4,[2;1]); RotationConstraint(link3,link4,ax)]
joint24 = SocketConstraint(link2,link4,[2;2])

### ZAC_BEGIN
## Closed Chain
links = [baseLink; link1; link2; link3; link4]
constraints = [jointb; jointb1; joint12; jointb3; joint34; joint24]

## Open Chain (2 double pendulums)
# links = [baseLink; link1; link2; link3; link4]
# constraints = [jointb; jointb1; joint12; jointb3; joint34]

## Single Pendulum
# links = [baseLink; link1]
# constraints = [jointb; jointb1]

## Single Body
# links = [baseLink]
# constraints = jointb

bot = Robot(links, constraints)
# bot = Robot(links)

simul = Simulation(bot)

# sim! is the main function. Call this with the debugger
# sim!(simul,debug=false,disp=true) # Set debug=true to see all warning messages
# trajS = simul.trajS

# include("visualizeTwoTwoBar.jl")  # To visualize the Closed or Open Chain with four links

### ZAC_END
