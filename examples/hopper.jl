using ConstrainedDynamics
using ConstrainedDynamicsVis

# Initial orientation
ϕ1 = 0;
q1 = Quaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
head = Sphere(0.1, 0.5, color = RGBA(1., 1., 0.))
leg1 = Cylinder(0.05, 0.8, 0.2, color = RGBA(0., 0., 0.))

fricsandineqs1 = Friction(leg1, [0;0;1.0], 0.2; p = [0;0;-0.4])
frics1 = fricsandineqs1[1]
ineqcs1 = fricsandineqs1[2]


joint0to1 = EqualityConstraint(Floating(origin, head))
joint1to2 = EqualityConstraint(Prismatic(head, leg1, [0;0;1.0]))

links = [head;leg1]
eqcs = [joint0to1;joint1to2]

mech = Mechanism(origin, links, eqcs, ineqcs1, [frics1], Δt = 0.001)

setPosition!(head,x = [0.;0;0.4])
setPosition!(head,leg1)

firstjump = [false]

function controller!(mechanism, k)
    x = minimalCoordinates(mechanism, joint1to2)[1]
    dx = minimalVelocities(mechanism, joint1to2)[1]
    q = minimalCoordinates(mechanism, joint0to1)[4:5]
    dq = minimalVelocities(mechanism, joint0to1)[4:5]

    if k < 300
        uleg1 = 300*(0-x)+10.0*(0-dx)
        uhead = zeros(6)
        if !firstjump[1] && leg1.state.xc[3]-0.4 <0.01
            uleg1-=125
        elseif !firstjump[1] && leg1.state.xc[3]-0.4 >= 0.01
            firstjump[1] = true
            uhead = 120.0*[0;0.0;0;1/sqrt(2);1/sqrt(2);0]
        end

        setForce!(mechanism, joint0to1, uhead)
        setForce!(mechanism, joint1to2, [uleg1])
    elseif k>800
        uleg1 = 200*(0.0-x)+10.0*(0-dx)
        uhead = [0;0;0; 10*(0-q[1])+0.5*(0-dq[1]);10*(0-q[2])+0.5*(0-dq[2]);0]

        setForce!(mechanism, joint0to1, uhead)
        setForce!(mechanism, joint1to2, [uleg1])
    end
end

steps = Base.OneTo(1100)
storagehopper = Storage(steps,length(mech.bodies))

storagehopper = simulate!(mech, storagehopper, controller!, record = true, ε=1e-6)
visualize(mech, storagehopper)

# plot([g(ineqcs1[1].constraints[1],storage.x[2][i],storage.q[2][i])[1] for i=1:1000])


# function init()
#     setPosition!(link1,x = [0.;0;1.0])
#     ωtemp = (rand(3) .- 0.5) * 10
#     setVelocity!(link1,v = [0;2;4.],ω = ωtemp)
# end


