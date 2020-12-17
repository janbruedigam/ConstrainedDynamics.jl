using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width, depth = 0.5, 0.5
box1 = Box(width, depth, length1, 1., color = RGBA(1., 1., 0.))
b1 = Box(0.1, 0.1, .1, .1, color = RGBA(1., 0., 0.))

corners = [
    [[length1 / 2;length1 / 2;-length1 / 2]]
    [[length1 / 2;-length1 / 2;-length1 / 2]]
    [[-length1 / 2;length1 / 2;-length1 / 2]]
    [[-length1 / 2;-length1 / 2;-length1 / 2]]
    [[length1 / 2;length1 / 2;length1 / 2]]
    [[length1 / 2;-length1 / 2;length1 / 2]]
    [[-length1 / 2;length1 / 2;length1 / 2]]
    [[-length1 / 2;-length1 / 2;length1 / 2]]
]

# Links
origin = Origin{Float64}()

link1 = Body(box1)
links = [link1]

# Constraints
ineqcs = [InequalityConstraint(Impact(link1, [0;0;1.0], p = corners[i])) for i = 1:8]

eqcs = [EqualityConstraint(Floating(origin, link1))]

shapes = [box1;b1]


mech = Mechanism(origin, links, eqcs, ineqcs, shapes = shapes)
setPosition!(link1,x = [0.;-2;1.5])

setVelocity!(link1,v = [0;3;7.],ω = [0.1;0.1;0.1])
