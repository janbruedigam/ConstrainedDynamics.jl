using ConstrainedDynamics
using ConstrainedDynamicsVis
using BenchmarkTools

# Parameters
ex = [1.;0.;0.]

h = 1.
r = .25

vert11 = [0.;0.;h]
vert12 = -vert11*0

# Initial orientation
phi = pi / 4
q1 = Quaternion(RotX(phi))

function reset!(origin, links, N)
    setPosition!(origin,links[1],Δq = Quaternion(RotX(pi/2+pi/4)))
    setVelocity!(links[1])
    for j=2:N
        setPosition!(links[j-1],links[j],p1 = vert12, p2 = vert11)
        setVelocity!(links[j])
    end
end

mindepth = ones(50)
timing = zeros(50)

for N=1:50
    display(N)

    origin = Origin{Float64}()
    links = [Sphere(r, h, color=RGBA(0,0.3,0.8)) for i = 1:N]

    # Constraints
    jointb1 = EqualityConstraint(Floating(origin, links[1]))
    if N>1
        constraints = [jointb1;[EqualityConstraint(Spherical(links[i - 1], links[i]; p1=vert12, p2=vert11)) for i = 2:N]]
    else
        constraints = [jointb1]
    end

    fricsandineqs = [Friction(links[i], [0;0;1.0], 1.0) for i=1:N]
    frics = getindex.(fricsandineqs,1)
    ineqcs = vcat(getindex.(fricsandineqs,2)...)

    mech = Mechanism(origin, links, constraints, ineqcs, frics)

    storage = Storage(500,N)

    function controller!(mechanism, k)
        if k>2
            for i=1:N
                mindepth[N] = minimum([links[i].state.xc[3]; mindepth[N]])
            end
        end
    end

    # reset!(origin, links, N)
    # simulate!(mech, storage, record=true, ε=1e-3)
    # visualize(mech, storage)

    bm = @benchmarkable simulate!($mech, $storage, record=false, ε=1e-3) setup=(reset!($origin, $links, $N))
    timing[N] = BenchmarkTools.minimum(run(bm,samples=100,seconds=100)).time
end