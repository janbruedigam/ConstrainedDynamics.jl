mutable struct Simulation{T}
    t::T
    tend::T
    dt::T

    trajS::Array{T,2}

    robot::Robot{T}
end

function Base.show(io::IO, simulation::Simulation{T}) where T
    heading = string(simulation.tend,"s Simulation{",T,"} with time step dt = ",simulation.dt,":")
    robot = string("\n Robot (robot): ",simulation.robot)

    print(io,heading,robot)
end

function Simulation(robot::Robot{T,Nl}; t=0, tend=10) where {T,Nl}
    nP = 3
    trajS = zeros(T,Nl*(nP+1),Int(ceil((tend-t)/robot.dt))+1)
    Simulation{T}(t,tend,robot.dt,trajS,robot)
end


function sim!(simulation::Simulation{T};debug::Bool=false,disp::Bool=false) where T
    robot = simulation.robot
    Nl = length(robot.links)
    nLl = robot.nLl
    nP = 3

    counter = 1;
    while simulation.t<simulation.tend
        x, φ = getState(robot)
        simulation.trajS[1:Nl*nP,counter] = x
        simulation.trajS[Nl*nP+1:Nl*(nP+1),counter] = φ

        ### ZAC_BEGIN
        sol = newton(robot,warning=debug)[1]
        updateRobot!(robot,sol)
        ### ZAC_END

        counter+=1;
        simulation.t += simulation.dt;
        if disp && (counter*simulation.dt)%1<simulation.dt*(1-.1)
            display(counter*simulation.dt)
        end
        # break
    end

    # return trajS
end


function plotTraj(simulation,trajS,id)
    p = plot(collect(0:simulation.dt:simulation.tend),trajS[id[1],:])
    for element in Iterators.rest(id,2)
        plot!(collect(0:simulation.dt:simulation.tend),trajS[element,:])
    end
    return p
end
