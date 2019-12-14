using Blink
using Colors: RGBA, RGB # Handle RGB colors
using CoordinateTransformations # Translations and rotations
using FileIO # Save and load files
using GeometryTypes: # Define geometric shapes
    GeometryTypes, HyperRectangle, Vec, Point, Rectangle, Cylinder,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh, Pyramid
using LinearAlgebra
using MeshCat # Visualize 3D animations
using MeshIO # Load meshes in MeshCat

vis = Visualizer()
open(vis, Blink.Window())

radius = 0.05
l1 = 1.
l2 = sqrt(2)/2
vislink1 = Cylinder(Point(0.0,0.0,-l1/2),Point(0.0,0.0,l1/2), radius)
vislink2 = Cylinder(Point(0.0,0.0,-l2/2),Point(0.0,0.0,l2/2), radius)
vislink3 = Cylinder(Point(0.0,0.0,-l1/2),Point(0.0,0.0,l1/2), radius)
vislink4 = Cylinder(Point(0.0,0.0,-l2/2),Point(0.0,0.0,l2/2), radius)

yellow_color = [1, 1, 0]
red_color = [1, 0, 0]
yellow_mat = MeshPhongMaterial(color=RGBA(yellow_color...))
red_mat = MeshPhongMaterial(color=RGBA(red_color...))

setobject!(vis["bundle/vislink1"], vislink1, red_mat)
setobject!(vis["bundle/vislink2"], vislink2, red_mat)
setobject!(vis["bundle/vislink3"], vislink3, yellow_mat)
setobject!(vis["bundle/vislink4"], vislink4, yellow_mat)

settransform!(vis["bundle/vislink1"], compose(Translation(bot.storage.x[1][1]...),LinearMap(Quat(bot.storage.q[1][1]...))))
settransform!(vis["bundle/vislink2"], compose(Translation(bot.storage.x[2][1]...),LinearMap(Quat(bot.storage.q[2][1]...))))
settransform!(vis["bundle/vislink3"], compose(Translation(bot.storage.x[3][1]...),LinearMap(Quat(bot.storage.q[3][1]...))))
settransform!(vis["bundle/vislink4"], compose(Translation(bot.storage.x[4][1]...),LinearMap(Quat(bot.storage.q[4][1]...))))

framerate = 100
anim = MeshCat.Animation(Dict{MeshCat.SceneTrees.Path,MeshCat.AnimationClip}(), framerate)

N = 1000 # Number of frames
for k=1:N
    # Set the pose of the car bundle.
    MeshCat.atframe(anim, vis, k) do frame
        settransform!(vis["bundle/vislink1"], compose(Translation(bot.storage.x[1][k]...),LinearMap(Quat(bot.storage.q[1][k]...))))
        settransform!(vis["bundle/vislink2"], compose(Translation(bot.storage.x[2][k]...),LinearMap(Quat(bot.storage.q[2][k]...))))
        settransform!(vis["bundle/vislink3"], compose(Translation(bot.storage.x[3][k]...),LinearMap(Quat(bot.storage.q[3][k]...))))
        settransform!(vis["bundle/vislink4"], compose(Translation(bot.storage.x[4][k]...),LinearMap(Quat(bot.storage.q[4][k]...))))
    end
end

MeshCat.setanimation!(vis, anim);
