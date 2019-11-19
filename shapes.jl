box(x,y,z,m::T) where T = m, 1/12*m*diagm([y^2+z^2;x^2+z^2;x^2+y^2])
