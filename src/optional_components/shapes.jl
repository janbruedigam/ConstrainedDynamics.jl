abstract type Shape{T} end

mutable struct Mesh{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    path::String
    color::RGBA

    function Mesh(path::String, m::Real, J::AbstractMatrix;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((m, J, xoffset, qoffset))...)

        bodyids = Int64[]
        new{T}(bodyids, m, J, xoffset, qoffset, path, color)
    end
end

mutable struct Box{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    xyz::SVector{3,T}
    color::RGBA

    function Box(x::Real, y::Real, z::Real, m::Real;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((x, y, z, m, xoffset, qoffset))...)

        J = 1 / 12 * m * diagm([y^2 + z^2;x^2 + z^2;x^2 + y^2])
        bodyids = Int64[]

        new{T}(bodyids, m, J, xoffset, qoffset, [x;y;z], color)
    end
end

mutable struct Cylinder{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    rh::SVector{2,T}
    color::RGBA

    # Cylinder points in the z direction
    function Cylinder(r::Real, h::Real, m::Real;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((r, h, m, xoffset, qoffset))...)

        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2;r^2 + 1 / 6 * h^2;r^2])
        bodyids = Int64[]

        new{T}(bodyids, m, J, xoffset, qoffset, [r;h], color)
    end
end

mutable struct Sphere{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    r::T
    color::RGBA

    function Sphere(r::Real, m::Real;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((r, m, xoffset, qoffset))...)
        
        J = 2 / 5 * m * diagm([r^2 for i = 1:3])
        bodyids = Int64[]

        new{T}(bodyids, m, J, xoffset, qoffset, r, color)
    end
end

