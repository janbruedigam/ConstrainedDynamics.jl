Spherical(body1,body2,p1,p2) = Translational0(body1,body2,p1,p2)
Cylindrical(body1,body2,p1,p2,axis) = Translational1(body1,body2,p1,p2,axis)
Planar(body1,body2,p1,p2,axis) = Translational2(body1,body2,p1,p2,axis)

Revolute(body1,body2,p1,p2,axis) = Translational0(body1,body2,p1,p2),Axis(body1,body2,axis)
