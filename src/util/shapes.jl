abstract type Shape{T} end

struct CustomShape{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xyz::SVector{3,T}
    color::RGBA

    # Currently only works for box shaped shapes
    # TODO use vectex description to create some kind of polygon
    # function CustomShape(m::T, J::AbstractMatrix{T}, vertices::Vector{AbstractVector{T}};color = RGBA(0.5, 0.5, 0.5)) where T
    function CustomShape(m::T, J::AbstractMatrix{T}, xyz::AbstractVector{T};color = RGBA(0.5, 0.5, 0.5)) where T
        bodyids = Int64[]
        new{T}(bodyids, m, J, xyz, color)
    end
end

struct Box{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xyz::SVector{3,T}
    color::RGBA

    function Box(x, y, z, m::T;color = RGBA(0.5, 0.5, 0.5)) where T
        J = 1 / 12 * m * diagm([y^2 + z^2;x^2 + z^2;x^2 + y^2])
        bodyids = Int64[]

        new{T}(bodyids, m, J, [x;y;z], color)
    end
end

struct Cylinder{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    rh::SVector{2,T}
    color::RGBA

    # by default the cylinder points in the z direction
    # TODO other direction
    function Cylinder(r, h, m::T;color = RGBA(0.5, 0.5, 0.5)) where T
        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2;r^2 + 1 / 6 * h^2;r^2])
        bodyids = Int64[]

        new{T}(bodyids, m, J, [r;h], color)
    end
end

