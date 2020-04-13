using MeshCat
using Blink
vis = Visualizer()
open(vis, Blink.Window())

using GeometryTypes
using CoordinateTransformations

setobject!(vis, HyperRectangle(Vec(0., 0, 0), Vec(1., 1, 1)))
settransform!(vis, Translation(-0.5, -0.5, 0.5))