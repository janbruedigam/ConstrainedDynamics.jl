# t3r3: fixed
Fixed(body1::AbstractBody{T},body2,p1,p2;offset = Quaternion{T}()) where T = Translational3{T}(body1, body2, p1=p1, p2=p2), Rotational3{T}(body1, body2, offset = offset)

# t2r3: prismatic
Prismatic(body1::AbstractBody{T},body2,p1,p2,axis;offset = Quaternion{T}()) where T = Translational2{T}(body1, body2, p1=p1, p2=p2, axis=axis), Rotational3{T}(body1, body2, offset = offset)

# t1r3: ?
# ?

# t0r3: fixed orientation (chicken's head)
FixedOrientation(body1::AbstractBody{T},body2;offset = Quaternion{T}()) where T = Translational0{T}(body1, body2), Rotational3{T}(body1, body2, offset = offset)

# t3r2: Revolute
Revolute(body1::AbstractBody{T},body2,p1,p2,axis;offset = Quaternion{T}()) where T = Translational3{T}(body1, body2, p1=p1, p2=p2), Rotational2{T}(body1, body2, axis = axis, offset = offset)

# t2r2: Cylindrical
Cylindrical(body1::AbstractBody{T},body2,p1,p2,axis;offset = Quaternion{T}()) where T = Translational2{T}(body1, body2, p1=p1, p2=p2, axis=axis), Rotational2{T}(body1, body2, axis = axis, offset = offset)

# t1r2: ?
# ?

# t0r2: ?
# ?

# What is r1?
# t3r1: universal?
# ?
# t2r1: ?
# ?
# t1r1: ?
# ?
# t0r1: ?
# ?

# t3r0: Spherical
Spherical(body1::AbstractBody{T},body2,p1,p2) where T = Translational3{T}(body1, body2, p1=p1, p2=p2), Rotational0{T}(body1, body2)

# t2r0: Point-on-Line
CylindricalFree(body1::AbstractBody{T},body2,p1,p2,axis) where T = Translational2{T}(body1, body2, p1=p1, p2=p2, axis=axis), Rotational0{T}(body1, body2)

# t1r0: Planar
Planar(body1::AbstractBody{T},body2,p1,p2,axis) where T = Translational1{T}(body1, body2, p1=p1, p2=p2, axis=axis), Rotational0{T}(body1, body2)

# t0r0: Free
OriginConnection(body1::AbstractBody{T},body2) where T = Translational0{T}(body1, body2), Rotational0{T}(body1, body2)
