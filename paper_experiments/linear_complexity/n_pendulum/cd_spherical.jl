using ConstrainedDynamics
using ConstrainedDynamicsVis
using BenchmarkTools


# Parameters
ex = [1.;0.;0.]

h = 1.
r = .05

vert11 = [0.;0.;h / 2]
vert12 = -vert11

# Initial orientation
phi = pi / 4
q1 = Quaternion(RotX(phi))

function reset!(mechanism, origin, links)
    setPosition!(origin,links[1],p2 = vert11,Δq = q1)
    previd = links[1].id
    for body in Iterators.drop(mechanism.bodies, 1)
        setPosition!(ConstrainedDynamics.getbody(mechanism, previd), body, p1 = vert12, p2 = vert11)
        previd = body.id
        setVelocity!(body)
    end
end

timing = zeros(100)

for N=1:100
    display(N)
    # Links

    origin = Origin{Float64}()
    links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:N]

    # Constraints
    jointb1 = EqualityConstraint(Spherical(origin, links[1]; p2=vert11))
    if N>1
        constraints = [jointb1;[EqualityConstraint(Spherical(links[i - 1], links[i]; p1=vert12, p2=vert11)) for i = 2:N]]
    else
        constraints = [jointb1]
    end


    mech = Mechanism(origin, links, constraints;Δt = 0.01)

    storage = Storage(1000,N)

    bm = @benchmarkable simulate!($mech,$storage,record=false) setup=(reset!($mech, $origin, $links))
    timing[N] = BenchmarkTools.minimum(run(bm,samples=100,seconds=100)).time

end
# storage = simulate!(mech, 10., record = true)
# visualize(mech, storage)