using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]

l1 = 1.0
x, y = .1, .1
b1 = Box(x, y, l1, l1, color = RGBA(1., 1., 0.))

vert11 = [0.;0.;l1 / 2]
vert12 = -vert11
verts = [[vert11];[vert12]]

# Initial orientation
offset1 = pi / 4
offset2 = pi / 2
phi1 = pi / 8
q1 = Quaternion(RotX(phi1))
qoff1 = Quaternion(RotX(offset1))
qoff2 = Quaternion(RotX(offset2))

N = 10

# Links
origin = Origin{Float64}()
links = [Body(b1) for i = 1:4 * N]


function fourbar(links, vertices, axis)
    j1 = EqualityConstraint(Revolute(links[1], links[2], vertices[1], vertices[2], axis))
    j2 = EqualityConstraint(Revolute(links[2], links[3], vertices[3], vertices[2], axis), Cylindrical(links[2], links[4], vertices[2], vertices[2], axis))
    j3 = EqualityConstraint(Revolute(links[4], links[5], vertices[3], vertices[2], axis))
    j4 = EqualityConstraint(Revolute(links[3], links[5], vertices[3], vertices[3], axis))

    return j1, j2, j3, j4
end

function initfourbar!(mechanism, links, vertices, Δq1, Δq2)
    setPosition!(mechanism, links[1], links[2], p1 = vertices[1], p2 = vertices[2], Δq = Δq1)
    setPosition!(mechanism, links[2], links[3], p1 = vertices[3], p2 = vertices[2], Δq = inv(Δq2) * inv(Δq2))
    setPosition!(mechanism, links[2], links[4], p1 = vertices[2], p2 = vertices[2], Δq = inv(Δq2) * inv(Δq2))
    setPosition!(mechanism, links[4], links[5], p1 = vertices[3], p2 = vertices[2], Δq = Δq2 * Δq2)
end

# Constraints
constraints = [fourbar([origin;links[1:4]], [[zeros(3)];verts], ex)...]
for i = 2:N
    push!(constraints, fourbar(links[(i - 1) * 4:i * 4], [[vert12];verts], ex)...)
end


shapes = [b1]

mech = Mechanism(origin, links, constraints, shapes = shapes)
if N > 1
    initfourbar!(mech, [origin;links[1:4]], [[zeros(3)];verts], q1 * qoff1, q1)
else
    initfourbar!(mech, [origin;links[1:4]], [[zeros(3)];verts], q1 * qoff2, q1)
end

for i = 2:N - 1
    initfourbar!(mech, links[(i - 1) * 4:i * 4], [[vert12];verts], q1 / q1, q1)
end

if N > 1
    initfourbar!(mech, links[(N - 1) * 4:N * 4], [[vert12];verts], inv(qoff1) * qoff2, q1)
end


simulate!(mech,save = true)
visualize!(mech)
