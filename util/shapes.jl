abstract type Shape end

mutable struct Box{T} <: Shape
    linkids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    lwh::SVector{3,T}
    color::RGBA

    function Box(x,y,z,m::T;color=RGBA(0.5,0.5,0.5)) where T
        J = 1/12*m*diagm([y^2+z^2;x^2+z^2;x^2+y^2])
        linkids = Vector{Int64}(undef,0)

        new{T}(linkids,m,J,[x;y;z],color)
    end
end
