include("dynamics.jl")
include("control.jl")


rng = MersenneTwister(1)
M = 100
paramcontainer = [[0.1; 0; 1; 0; -1.5]]
paramstorage = [[0.1; 0; 1; 0; -1.5]]
bias = zeros(5)
distance = 0.0
explore_factor = 0.1
distancestorage = zeros(M)

function reset_state!()
    setPosition!(mech, floating_constraint, zeros(6))
    for i=1:4
        setPosition!(mech, legjoints[(i-1)*3+1],[0])
        setPosition!(mech, legjoints[(i-1)*3+2],[paramcontainer[1][3]])
        setPosition!(mech, legjoints[(i-1)*3+3],[paramcontainer[1][5]])
    end
    setPosition!(mech, floating_constraint,[0, 0, -calf_body.state.xc[3]-ConstrainedDynamics.vrotate(contact_point,calf_body.state.qc)[3] + 0.01, 0, 0, 0])

    setVelocity!(mech, floating_constraint, zeros(6))
    for i=1:12
        setVelocity!(mech, legjoints[i],[0])
    end
end

for i=1:M
    display("run: $i")

    if bias == zeros(5)
        paramcontainer[1] += randn!(rng, zeros(5))*explore_factor
    else
        paramcontainer[1] += randn!(rng, zeros(5))*0.002 + normalize(bias)*0.01
    end

    reset_state!()
    x0 = trunk_body.state.xc[1]

    try
        storagequad = simulate!(mech, storagequad, controller!, ε = 1e-3)
    catch
        display("  errored")
        paramcontainer[1] = paramstorage[end]
        bias = zeros(5)
        explore_factor *= 0.9
    end

    distancenew = trunk_body.state.xc[1] - x0

    if any(isnan.(getindex.(storagequad.x[1],1))) || distancenew <= distance
        display("  unsuccessful")
        any(isnan.(getindex.(storagequad.x[1],1))) && display("  nans")
        paramcontainer[1] = paramstorage[end]
        bias = zeros(5)
        explore_factor *= 0.9
    else
        display("  successful")
        distance = distancenew
        push!(paramstorage,paramcontainer[1])
        bias = paramstorage[end]-paramstorage[end-1]
        explore_factor = 0.1
    end

    display("  distance: $distancenew")
    distancestorage[i] = distance
end

reset_state!()

storagequad = simulate!(mech, storagequad, controller!, ε = 1e-3)
visualize(mech, storagequad, showframes = false)