abstract type GenericJoint{T,N} <: AbstractJoint{T,N} end

@inline constraintmat(::GenericJoint) = I

g(joint::GenericJoint, statea::State, stateb::State, Δt) = joint.g(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
g(joint::GenericJoint, stateb::State, Δt) = joint.g(joint, posargsnext(stateb, Δt)...)
g(joint::GenericJoint, statea::State, stateb::State) = joint.g(joint, posargsc(statea)..., posargsc(stateb)...)
g(joint::GenericJoint, stateb::State) = joint.g(joint, posargsc(stateb)...)