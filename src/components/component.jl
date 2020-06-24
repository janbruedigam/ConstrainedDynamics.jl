abstract type Component{T} end
abstract type AbstractBody{T} <: Component{T} end
abstract type AbstractConstraint{T,N} <: Component{T} end

# TODO do id differently?
CURRENTID = -1
getGlobalID() = (global CURRENTID -= 1; return CURRENTID + 1)
resetGlobalID() = (global CURRENTID = -1; return)

Base.show(io::IO, component::Component) = summary(io, component)

Base.eltype(::Type{<:Component{E}}) where {E} = @isdefined(E) ? E : Any
Base.length(::AbstractConstraint{T,N}) where {T,N} = N
getid(component::Component) = component.id
activate!(component::Component) = (component.active = true; return)
deactivate!(component::Component) = (component.active = false; return)

isactive(component::Component) = component.active
isinactive(component::Component) = !isactive(component)
