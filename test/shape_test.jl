using ConstrainedDynamics
using Rotations

box = Box(rand(4)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(UnitQuaternion))
Body(box)
@test true
cylinder = Cylinder(rand(3)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(UnitQuaternion))
Body(cylinder)
@test true
sphere = Sphere(rand(2)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(UnitQuaternion))
Body(sphere)
@test true
mesh = Mesh("test.obj", rand(), rand(3,3); color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(UnitQuaternion))
Body(mesh)
@test true