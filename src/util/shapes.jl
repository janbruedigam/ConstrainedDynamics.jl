abstract type Shape{T} end

mutable struct Mesh{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoff::SVector{3,T}
    qoff::UnitQuaternion{T}

    path::String
    color::RGBA

    function Mesh(path::String, m::Real, J::AbstractMatrix;
            color = RGBA(0.75, 0.75, 0.75), xoff::AbstractVector = zeros(3), qoff::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((m, J, xoff, qoff))...)

        bodyids = Int64[]
        new{T}(bodyids, m, J, xoff, qoff, path, color)
    end
end

mutable struct Box{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoff::SVector{3,T}
    qoff::UnitQuaternion{T}

    xyz::SVector{3,T}
    color::RGBA

    function Box(x::Real, y::Real, z::Real, m::Real;
            color = RGBA(0.75, 0.75, 0.75), xoff::AbstractVector = zeros(3), qoff::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((x, y, z, m, xoff, qoff))...)

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
    qoff::UnitQuaternion{T}

    rh::SVector{2,T}
    color::RGBA

    # Cylinder points in the z direction
    function Cylinder(r::Real, h::Real, m::Real;
            color = RGBA(0.75, 0.75, 0.75), xoff::AbstractVector = zeros(3), qoff::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((r, h, m, xoff, qoff))...)

        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2;r^2 + 1 / 6 * h^2;r^2])
        bodyids = Int64[]

        new{T}(bodyids, m, J, xoff, qoff, [r;h], color)
    end
end

mutable struct Sphere{T} <: Shape{T}
    bodyids::Vector{Int64}
    m::T
    J::SMatrix{3,3,T,9}

    xoff::SVector{3,T}
    qoff::UnitQuaternion{T}

    r::T
    color::RGBA

    function Sphere(r::Real, m::Real;
            color = RGBA(0.75, 0.75, 0.75), xoff::AbstractVector = zeros(3), qoff::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((r, h, m, xoff, qoff))...)
        
        J = 2 / 5 * m * diagm([r^2 for i = 1:3])
        bodyids = Int64[]

        new{T}(bodyids, m, J, xoff, qoff, r, color)
    end
end

