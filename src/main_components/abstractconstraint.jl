@inline function ∂g∂ʳpos(mechanism, eqc::AbstractConstraint, id::Integer)
    id == eqc.parentid ? (return ∂g∂ʳposa(mechanism, eqc, id)) : (return ∂g∂ʳposb(mechanism, eqc, id))
end
@inline function ∂g∂ʳvel(mechanism, eqc::AbstractConstraint, id::Integer)
    id == eqc.parentid ? (return ∂g∂ʳvela(mechanism, eqc, id)) : (return ∂g∂ʳvelb(mechanism, eqc, id))
end