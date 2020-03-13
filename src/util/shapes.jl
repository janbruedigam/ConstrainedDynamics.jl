abstract type Shape{T} end

mutable struct Mesh{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoff::SVector{3,T}
    qoff::Quaternion{T}

    path::String
    color::RGBA

    function Mesh(path::String, m::T, J::AbstractMatrix{T};color = RGBA(0.75, 0.75, 0.75), xoff::AbstractVector{T}=zeros(3), qoff::Quaternion{T}=Quaternion{T}()) where T
        bodyids = Int64[]
        new{T}(bodyids, m, J, xoff, qoff, path, color)
    end
end

mutable struct Box{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoff::SVector{3,T}
    qoff::Quaternion{T}

    xyz::SVector{3,T}
    color::RGBA

    function Box(x, y, z, m::T;color = RGBA(0.75, 0.75, 0.75), xoff::AbstractVector{T}=zeros(3), qoff::Quaternion{T}=Quaternion{T}()) where T
        J = 1 / 12 * m * diagm([y^2 + z^2;x^2 + z^2;x^2 + y^2])
        bodyids = Int64[]

        new{T}(bodyids, m, J, xoff, qoff, [x;y;z], color)
    end
end

mutable struct Cylinder{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoff::SVector{3,T}
    qoff::Quaternion{T}

    rh::SVector{2,T}
    color::RGBA

    # by default the cylinder points in the z direction
    # TODO other direction
    function Cylinder(r, h, m::T;color = RGBA(0.75, 0.75, 0.75), xoff::AbstractVector{T}=zeros(3), qoff::Quaternion{T}=Quaternion{T}()) where T
        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2;r^2 + 1 / 6 * h^2;r^2])
        bodyids = Int64[]

        new{T}(bodyids, m, J, xoff, qoff, [r;h], color)
    end
end

