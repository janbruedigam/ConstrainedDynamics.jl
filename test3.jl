using Rotations
using BenchmarkTools
using TimerOutputs
using Plots

(@isdefined FullCordDynamics) ? nothing : include("FullCordDynamics.jl")
using Main.FullCordDynamics

bot = parse_urdf("twoTwoBarDiffLength.urdf")
origin = bot.root
links = bot.nodes[bot.nodesRange[1]]
constraints = bot.nodes[bot.nodesRange[2]]

phi1, phi2, phi3, phi4 = pi/2, -pi/4, 0., 3*pi/4.
q1, q2, q3, q4 = Quaternion(RotX(phi1)), Quaternion(RotX(phi2)), Quaternion(RotX(phi3)), Quaternion(RotX(phi4))

setInit!(origin,links[1],[2;1],q=q1)
setInit!(links[1],links[2],[2;1],q=q2)
setInit!(origin,links[3],[3;1],q=q3)
setInit!(links[3],links[4],[2;1],q=q4)

sim!(bot,save=true)
trajS = trajSFunc(bot)
# include("visualizeTwoTwoBar.jl")
