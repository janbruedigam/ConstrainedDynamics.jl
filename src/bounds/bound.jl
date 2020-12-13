abstract type Bound{T} end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, bound::Bound{T}) where {T}
    summary(io, bound)
end

### General functions
@inline getT(bound::Bound{T}) where T = T


### Constraints and derivatives
## Position level constraint wrappers
g(contact::Bound, body::Body, Δt) = g(contact, body.state, Δt)

## Additional force for friction default
@inline additionalforce(bound::Bound{T}) where T = szeros(T, 6)


# @inline g(bound::Bound{T}) where T = zero(T)

# @inline ∂g∂pos(bound::Bound{T}) where T = szeros(T, 6)
# @inline ∂g∂vel(bound::Bound{T}) where T = szeros(T, 6)
# @inline schurf(bound::Bound{T}) where T = szeros(T, 6)
# @inline schurD(bound::Bound{T}) where T = szeros(T, 6, 6)



