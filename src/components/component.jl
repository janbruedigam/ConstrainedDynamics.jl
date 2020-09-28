abstract type Component{T} end
abstract type AbstractBody{T} <: Component{T} end
abstract type AbstractConstraint{T,N} <: Component{T} end
abstract type AbstractJoint{T,N} end

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

function Base.getindex(dict::UnitDict{Base.OneTo{K},<:Component}, key::String) where K
    for component in dict.values
        component.name == key && return component
    end
    
    return
end
function Base.getindex(dict::UnitDict{UnitRange{K},<:Component}, key::String) where K
    for component in dict.values
        component.name == key && return component
    end
    
    return
end