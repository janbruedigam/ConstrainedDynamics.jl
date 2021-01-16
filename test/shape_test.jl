using ConstrainedDynamics
using Rotations

Box(rand(4)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(UnitQuaternion))
@test true
Cylinder(rand(3)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(UnitQuaternion))
@test true
Sphere(rand(2)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(UnitQuaternion))
@test true
Mesh("test.obj", rand(), rand(3,3); color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(UnitQuaternion))
@test true