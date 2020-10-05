abstract type GenericJoint{T,N} <: AbstractJoint{T,N} end

@inline g(joint::GenericJoint, body1::Body, body2::Body, Δt) =  g(joint, body1.state, body2.state, Δt) 
@inline g(joint::GenericJoint, body1::Origin, body2::Body, Δt) = g(joint, body2.state, Δt)            
@inline g(joint::GenericJoint, body1::Body, body2::Body) =  g(joint, body1.state, body2.state)               
@inline g(joint::GenericJoint, body1::Origin, body2::Body) = g(joint, body2.state)                          

@inline function ∂g∂ʳposa(joint::GenericJoint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return  ∂g∂ʳposa(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposa(joint)
    end
end
@inline function ∂g∂ʳposb(joint::GenericJoint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return  ∂g∂ʳposb(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end
@inline function ∂g∂ʳposb(joint::GenericJoint, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return  ∂g∂ʳposb(joint, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end

@inline function ∂g∂posac(joint::GenericJoint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posac(joint, body1.state, body2.state)
    else
        return ∂g∂posac(joint)
    end
end
@inline function ∂g∂posbc(joint::GenericJoint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return  ∂g∂posbc(joint, body1.state, body2.state)
    else
        return ∂g∂posbc(joint)
    end
end
@inline function ∂g∂posbc(joint::GenericJoint, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posbc(joint, body2.state)
    else
        return ∂g∂posbc(joint)
    end
end

@inline function ∂g∂ʳvela(joint::GenericJoint, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return  ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvela(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::GenericJoint, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return  ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::GenericJoint, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end

@inline function ∂Fτ∂ua(joint::GenericJoint, body1::Body)
    #return ∂Fτ∂ua(joint, body1.state) * zerodimstaticadjoint(nullspacemat(joint))
end
@inline function ∂Fτ∂ub(joint::GenericJoint, body1::Body, body2::Body)
    #=if body2.id == joint.childid
        return ∂Fτ∂ub(joint, body1.state, body2.state) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return ∂Fτ∂ub(joint)
    end=#
end
@inline function ∂Fτ∂ub(joint::GenericJoint, body1::Origin, body2::Body)
   #= if body2.id == joint.childid
        return return ∂Fτ∂ub(joint, body2.state) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return ∂Fτ∂ub(joint)
    end=#
end
g(joint::GenericJoint, statea::State, stateb::State, Δt) = joint.g(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
g(joint::GenericJoint, stateb::State, Δt) = joint.g(joint, posargsnext(stateb, Δt)...)
g(joint::GenericJoint, statea::State, stateb::State) = joint.g(joint, posargsc(statea)..., posargsc(stateb)...)
g(joint::GenericJoint, stateb::State) = joint.g(joint, posargsc(stateb)...)