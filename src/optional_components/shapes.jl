abstract type Shape{T} end

struct EmptyShape{T} <: Shape{T}
    EmptyShape() = new{Float64}()
end

mutable struct Mesh{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    path::String
    scale::SVector{3,T}
    color::RGBA

    function Mesh(path::String;
            scale::AbstractVector = ones(3), color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((xoffset, qoffset))...)

        new{T}(xoffset, qoffset, path, scale, color)
    end

    function Mesh(path::String, m::Real, J::AbstractMatrix;
            scale::AbstractVector = ones(3), name::String="", color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((m, J, xoffset, qoffset))...)

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, path, scale, color))
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
        T = promote_type(eltype.((x, y, z, xoffset, qoffset))...)

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
        T = promote_type(eltype.((r, h, xoffset, qoffset))...)

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
        T = promote_type(eltype.((r, xoffset, qoffset))...)

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

mutable struct Pyramid{T} <: Shape{T}
    xoffset::SVector{3,T}
    qoffset::UnitQuaternion{T}

    wh::SVector{2,T}
    color::RGBA

    # Pyramid points in the z direction, Center of mass at 1/4 h
    function Pyramid(w::Real, h::Real;
            color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((w, h, xoffset, qoffset))...)

        new{T}(xoffset, qoffset, [w;h], color)
    end

    function Pyramid(w::Real, h::Real, m::Real;
            name::String="", color = RGBA(0.75, 0.75, 0.75), xoffset::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion)
        )
        T = promote_type(eltype.((w, h, m, xoffset, qoffset))...)
        J = 1/80 * m * diagm([4*w^2+3*h^2;4*w^2+3*h^2;8*w^2])

        return Body(m, J; name=name, shape=new{T}(xoffset, qoffset, [w;h], color))
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
