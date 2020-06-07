using ConstrainedDynamics
using Rotations

box = Box(rand(4)...; color = RGBA(rand(3)...), xoff = rand(3), qoff = Quaternion(rand(RotMatrix{3})))
@test true
cylinder = Cylinder(rand(3)...; color = RGBA(rand(3)...), xoff = rand(3), qoff = Quaternion(rand(RotMatrix{3})))
@test true
sphere = Sphere(rand(2)...; color = RGBA(rand(3)...), xoff = rand(3), qoff = Quaternion(rand(RotMatrix{3})))
@test true
mesh = Mesh("test.obj", rand(), rand(3,3); color = RGBA(rand(3)...), xoff = rand(3), qoff = Quaternion(rand(RotMatrix{3})))
@test true