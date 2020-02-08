# t0r0: fixed
Fixed(body1::AbstractBody{T},body2,p1,p2;offset=Quaternion{T}()) where T = Translational0(body1,body2,p1,p2),Rotational0(body1,body2,offset=offset)

# t1r0: prismatic
Prismatic(body1,body2,p1,p2,axis;offset=Quaternion{T}()) where T = Translational1(body1,body2,p1,p2,axis),Rotational0(body1,body2,offset=offset)

# t2r0: ?
# ?

# t3r0: fixed orientation (chicken's head)
FixedOrientation(body1::AbstractBody{T},body2;offset=Quaternion{T}()) where T = Rotational0(body1,body2,offset=offset)

# t0r1: Revolute
Revolute(body1,body2,p1,p2,axis) = Translational0(body1,body2,p1,p2),Rotational1(body1,body2,axis)

# t1r1: Cylindrical
Cylindrical(body1,body2,p1,p2,axis) = Translational1(body1,body2,p1,p2,axis),Rotational1(body1,body2,axis)

# t2r1: ?
# ?

# t3r1: ?
# ?

# t0r2: universal?
# ?

# t1r2: ?
# ?

# t2r2: ?
# ?

# t3r2: ?
# ?

# t0r3: Spherical
Spherical(body1,body2,p1,p2) = Translational0(body1,body2,p1,p2)

# t1r3: Point-on-Line
CylindricalFree(body1,body2,p1,p2,axis) = Translational1(body1,body2,p1,p2,axis)

# t2r3: Planar
Planar(body1,body2,p1,p2,axis) = Translational2(body1,body2,p1,p2,axis)

# t3r3: Free
OriginConnection(body1,body2) = TranslationalRotational6(body1,body2)
