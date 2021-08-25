@inline function ∂g∂ʳpos(mechanism, eqc::AbstractConstraint, body::Body)
    body.id == eqc.parentid ? (return ∂g∂ʳposa(mechanism, eqc, body)) : (return ∂g∂ʳposb(mechanism, eqc, body))
end
@inline function ∂g∂ʳvel(mechanism, eqc::AbstractConstraint, body::Body)
    body.id == eqc.parentid ? (return ∂g∂ʳvela(mechanism, eqc, body)) : (return ∂g∂ʳvelb(mechanism, eqc, body))
end