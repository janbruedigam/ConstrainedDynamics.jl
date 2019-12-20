abstract type Shape end

mutable struct Box{T} <: Shape
    linkids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xyz::SVector{3,T}
    color::RGBA

    function Box(x,y,z,m::T;color=RGBA(0.5,0.5,0.5)) where T
        J = 1/12*m*diagm([y^2+z^2;x^2+z^2;x^2+y^2])
        linkids = Vector{Int64}(undef,0)

        new{T}(linkids,m,J,[x;y;z],color)
    end
end

mutable struct Cylinder{T} <: Shape
    linkids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    rh::SVector{2,T}
    color::RGBA

    # by defaul the cylinder points in the z direction
    # TODO other direction
    function Cylinder(r,h,m::T;color=RGBA(0.5,0.5,0.5)) where T
        J = 1/2*m*diagm([r^2+1/6*h^2;r^2+1/6*h^2;r^2])
        linkids = Vector{Int64}(undef,0)

        new{T}(linkids,m,J,[r;h],color)
    end
end
