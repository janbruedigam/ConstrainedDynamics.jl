abstract type Component{T} end
abstract type AbstractBody{T} <: Component{T} end
abstract type AbstractConstraint{T,N} <: Component{T} end

# TODO do id differently?
CURRENTID = -1
getGlobalID() = (global CURRENTID -= 1; return CURRENTID + 1)
resetGlobalID() = (global CURRENTID = -1; return)

Base.show(io::IO, component::Component) = summary(io, component)

@inline Base.foreach(f, itr::Vector{<:Component}, arg...) = (for x in itr; f(x, arg...); end; return)

Base.length(::AbstractConstraint{T,N}) where {T,N} = N
activate!(constraint::AbstractConstraint) = (constraint.active = true; return)
deactivate!(constraint::AbstractConstraint) = (constraint.active = false; return)