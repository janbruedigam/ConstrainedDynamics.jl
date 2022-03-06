using ConstrainedDynamics
using ConstrainedDynamicsVis
using BenchmarkTools

# Initial orientation
ϕ1 = 0;
q1 = QuatRotation(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Cylinder(0.5,0.1, 1., color = RGBA(1., 0., 0.,0.5))

# data = zeros(50)
steps = Base.OneTo(100)

for i=15:15
display(i)

# Constraints
N = i

angle = [i*pi/N for i=1:N]
corners1 = [0.5*[sin(angle[i]);cos(angle[i]);0] for i=1:N]
corners2 = [0.5*[sin(angle[i]+pi);cos(angle[i]+pi);0] for i=1:N]
corners = [corners1;corners2]

fricsandineqs = [Friction(link1, [0;0;1.0], 1.0; p = corners[i]) for i=1:2*N]
frics = getindex.(fricsandineqs,1)
ineqcs = vcat(getindex.(fricsandineqs,2)...)


joint0to1 = EqualityConstraint(Floating(origin, link1))

links = [link1]
eqcs = [joint0to1]

mech = Mechanism(origin, links, eqcs, ineqcs, frics)

setPosition!(link1,x = [2.0;2.0;0.5],q=QuatRotation(RotY(pi/2-0.1)))
setVelocity!(link1,v=[0;-1*pi;0.0],ω = [0;0.0;2*pi])
storagedisc = simulate!(mech, 5, record = true, debug=false, ε=1e-6)
visualize(mech, storagedisc, showframes=true)

# storagedisc = Storage(steps,length(mech.bodies))

# t = @benchmarkable simulate!($mech, $steps, $storage, ε=1e-6)
# data[i] = BenchmarkTools.minimum(run(t,samples=100,seconds=100)).time/1e9
end