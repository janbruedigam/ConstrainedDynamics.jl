using ConstrainedDynamics
using Rotations: QuatRotation

# Body
Box(rand(4)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
Cylinder(rand(3)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
Sphere(rand(2)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
Pyramid(rand(3)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
Mesh("test.obj", rand(), rand(3,3); color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
# Just Shape
Box(rand(3)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
Cylinder(rand(2)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
Sphere(rand()...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
Pyramid(rand(2)...; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true
Mesh("test.obj"; color = RGBA(rand(3)...), xoffset = rand(3), qoffset = rand(QuatRotation).q)
@test true

storage1 = Storage{Float64}(Base.OneTo(2),2)
@test true
storage2 = Storage([[[[0.0;0;0] for i=1:2]];[[[0.0;0;0] for i=1:2]]],[[[Quaternion(RotX(0.0)) for i=1:2]];[[Quaternion(RotX(0.0)) for i=1:2]]])
@test true
