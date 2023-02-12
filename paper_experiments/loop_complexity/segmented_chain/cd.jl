using ConstrainedDynamics
using ConstrainedDynamicsVis


# Parameters
ex = [1.;0.;0.]

l1 = 1.0
x, y = .1, .1

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

function fourbar(links, vertices, axis)
    j1 = EqualityConstraint(Revolute(links[1], links[2], axis; p1=vertices[1], p2=vertices[2]))
    j2 = EqualityConstraint(Revolute(links[2], links[3], axis; p1=vertices[3], p2=vertices[2]), Cylindrical(links[2], links[4], axis; p1=vertices[2], p2=vertices[2]))
    j3 = EqualityConstraint(Revolute(links[4], links[5], axis; p1=vertices[3], p2=vertices[2]))
    j4 = EqualityConstraint(Revolute(links[3], links[5], axis; p1=vertices[3], p2=vertices[3]))

    return j1, j2, j3, j4
end

function initfourbar!(mechanism, links, vertices, Δq1, Δq2)
    setPosition!(links[1], links[2], p1 = vertices[1], p2 = vertices[2], Δq = Δq1)
    setPosition!(links[2], links[3], p1 = vertices[3], p2 = vertices[2], Δq = inv(Δq2) * inv(Δq2))
    setPosition!(links[2], links[4], p1 = vertices[2], p2 = vertices[2], Δq = inv(Δq2) * inv(Δq2))
    setPosition!(links[4], links[5], p1 = vertices[3], p2 = vertices[2], Δq = Δq2 * Δq2)
    links[1] != mechanism.origin && (links[1].state.vc = zeros(3))
    links[2].state.vc = zeros(3)
    links[3].state.vc = zeros(3)
    links[4].state.vc = zeros(3)
    links[5].state.vc = zeros(3)
    links[1] != mechanism.origin && (links[1].state.ωc = zeros(3))
    links[2].state.ωc = zeros(3)
    links[3].state.ωc = zeros(3)
    links[4].state.ωc = zeros(3)
    links[5].state.ωc = zeros(3)

    return
end

function reset!(mechanism, origin, links, N)
    if N > 1
        initfourbar!(mechanism, [origin;links[1:4]], [[zeros(3)];verts], q1 * qoff1, q1)
    else
        initfourbar!(mechanism, [origin;links[1:4]], [[zeros(3)];verts], q1 * qoff2, q1)
    end

    for i = 2:N - 1
        initfourbar!(mechanism, links[(i - 1) * 4:i * 4], [[vert12];verts], q1 / q1, q1)
    end

    if N > 1
        initfourbar!(mechanism, links[(N - 1) * 4:N * 4], [[vert12];verts], inv(qoff1) * qoff2, q1)
    end
end

timing = zeros(10)

for N = 1:10
    display(N)

    # Links
    origin = Origin{Float64}()
    links = [Box(x, y, l1, l1, color=RGBA(0,0.3,0.8)) for i = 1:4 * N]

    # Constraints
    constraints = [fourbar([origin;links[1:4]], [[zeros(3)];verts], ex)...]
    for i = 2:N
        push!(constraints, fourbar(links[(i - 1) * 4:i * 4], [[vert12];verts], ex)...)
    end


    mech = Mechanism(origin, links, constraints)

    storage = Storage(1000,N)

    # reset!(mech, origin, links, N)
    # storage = simulate!(mech, 10., record = true)
    # visualize(mech, storage)

    bm = @benchmarkable simulate!($mech,$storage,record=false) setup=(reset!($mech, $origin, $links, $N))
    timing[N] = BenchmarkTools.minimum(run(bm,samples=100,seconds=100)).time
end