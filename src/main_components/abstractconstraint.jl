@inline function ∂g∂ʳpos(mechanism, constraint::AbstractConstraint, body::Body)
    body.id == constraint.parentid ? (return ∂g∂ʳposa(mechanism, constraint, body)) : (return ∂g∂ʳposb(mechanism, constraint, body))
end
@inline function ∂g∂ʳvel(mechanism, constraint::AbstractConstraint, body::Body)
    body.id == constraint.parentid ? (return ∂g∂ʳvela(mechanism, constraint, body)) : (return ∂g∂ʳvelb(mechanism, constraint, body))
end