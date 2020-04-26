using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]

l1 = 1.0
l2 = 1.0# sqrt(2)/2
x, y = .1, .1
b1 = Box(x, y, l1, l1, color = RGBA(1., 1., 0.))
b2 = Box(x, y, l2, l2, color = RGBA(1., 1., 0.))
b3 = Box(x, y, l1, l1, color = RGBA(1., 0., 0.))
b4 = Box(x, y, l2, l2, color = RGBA(1., 0., 0.))

vert11 = [0.;0.;l1 / 2]
vert12 = -vert11

vert21 = [0.;0.;l2 / 2]
vert22 = -vert21

verts = [
    [vert11]
    [vert12]
    [vert21]
    [vert22]
]

# Initial orientation
phi1 = pi / 5
q1 = Quaternion(RotX(phi1))

# Links
origin = Origin{Float64}()

link1 = Body(b1)
link2 = Body(b2)
link3 = Body(b3)
link4 = Body(b4)

link5 = Body(b1)
link6 = Body(b2)
link7 = Body(b3)
link8 = Body(b4)

# Constraints
function fourbar(links, vertices, axis)
    j1 = EqualityConstraint(Revolute(links[1], links[2], zeros(3), vertices[1], axis))
    j2 = EqualityConstraint(Revolute(links[2], links[3], vertices[2], vertices[3], axis), Cylindrical(links[2], links[4], vertices[1], vertices[1], axis))
    j3 = EqualityConstraint(Revolute(links[4], links[5], vertices[2], vertices[3], axis))
    j4 = EqualityConstraint(Revolute(links[3], links[5], vertices[4], vertices[4], axis))

    return j1, j2, j3, j4
end

function initfourbar!(mechanism, links, vertices, Δq1, Δq2)
    setPosition!(mechanism, links[1], links[2], p2 = vertices[1], Δq = Δq1)
    setPosition!(mechanism, links[2], links[3], p1 = vertices[2], p2 = vertices[3], Δq = inv(Δq2))
    setPosition!(mechanism, links[2], links[4], p1 = vertices[1], p2 = vertices[1], Δq = inv(Δq2))
    setPosition!(mechanism, links[4], links[5], p1 = vertices[2], p2 = vertices[3], Δq = Δq2)
end

links = [link1; link2; link3; link4; link5; link6; link7; link8]
constraints = [fourbar([origin;links[1:4]], verts, ex)...,fourbar([links[4];links[5:8]], verts, ex)...]
shapes = [b1,b2,b3,b4]

mech = Mechanism(origin, links, constraints, InequalityConstraint{Float64}[], shapes = shapes)

initfourbar!(mech,[origin;links[1:4]],verts,q1,q1)
initfourbar!(mech,[links[4];links[5:8]],verts,q1 / q1,q1)


simulate!(mech,save = true)
visualize!(mech)
