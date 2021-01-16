abstract type Shape{T} end

struct EmptyShape{T} <: Shape{T}
    EmptyShape() = new{Float64}()
end

mutable struct Mesh{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    path::String
    color::RGBA


    function Mesh(path::String;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((xoffset, qoffset))...)

        new{T}(xoffset, qoffset, path, color)
    end

    function Mesh(path::String, m::Real, J::AbstractMatrix;
            name::String="", color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((m, J, xoffset, qoffset))...)

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, path, color))
    end
end

mutable struct Box{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    xyz::SVector{3,T}
    color::RGBA


    function Box(x::Real, y::Real, z::Real;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((x, y, z, m, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, [x;y;z], color)
    end

    function Box(x::Real, y::Real, z::Real, m::Real;
            name::String="", color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((x, y, z, m, xoffset, qoffset))...)
        J = 1 / 12 * m * diagm([y^2 + z^2;x^2 + z^2;x^2 + y^2])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, [x;y;z], color))
    end
end

mutable struct Cylinder{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    rh::SVector{2,T}
    color::RGBA

    # Cylinder points in the z direction
    function Cylinder(r::Real, h::Real;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((r, h, m, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, [r;h], color)
    end

    function Cylinder(r::Real, h::Real, m::Real;
            name::String="", color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((r, h, m, xoffset, qoffset))...)
        J = 1 / 2 * m * diagm([r^2 + 1 / 6 * h^2;r^2 + 1 / 6 * h^2;r^2])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, [r;h], color))
    end
end

mutable struct Sphere{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    r::T
    color::RGBA

    function Sphere(r::Real;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((r, m, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, r, color)
    end

    function Sphere(r::Real, m::Real;
            name::String="", color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((r, m, xoffset, qoffset))...)
        J = 2 / 5 * m * diagm([r^2 for i = 1:3])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, r, color))
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, shape::Shape{T}) where {T}
    summary(io, shape)
    println(io,"")
    println(io," xoffset: "*string(shape.xoffset))
    println(io," qoffset: "*string(shape.qoffset))
    println(io," color:   "*string(shape.color))
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, shape::EmptyShape{T}) where {T}
    summary(io, shape)
    println(io,"")
end
