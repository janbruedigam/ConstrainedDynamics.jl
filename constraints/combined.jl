include("../node.jl")

struct Combined{T,Nc,N,Nc²,NcN,Nl,C1,C2} <: Constraint{T,Nc,N,Nc²,NcN,Nl}
    constr1::C1
    constr2::C2

    data::NodeData{T,Nc,N,Nc²,NcN}
end

function Combined(constr1::C1, constr2::C2) where {C1,C2}
    # @assert constr1.T == constr2.T
    @assert linkids(constr1) == linkids(constr2)

    T = constr1.T
    Nc = constr1.Nc+constr2.Nc
    N = 6
    Nc² = Nc^2
    NcN = Nc*N
    Nl = length(linkids(constr1))
    data = NodeData{T,Nc,N}()

    Combined{T,Nc,N,Nc²,NcN,Nl,C1,C2}(constr1,constr2,data)
end


@inline g(C::Combined) = [g(C.constr1);g(C.constr2)]

@inline ∂g∂posa(C::Combined) = [∂g∂posa(C.constr1);∂g∂posa(C.constr2)]
@inline ∂g∂posb(C::Combined) = [∂g∂posb(C.constr1);∂g∂posb(C.constr2)]
@inline ∂g∂vela(C::Combined) = [∂g∂vela(C.constr1);∂g∂vela(C.constr2)]
@inline ∂g∂velb(C::Combined) = [∂g∂velb(C.constr1);∂g∂velb(C.constr2)]

@inline linkids(C::Combined{T,Nc,N,Nc²,NcN,1}) where {T,Nc,N,Nc²,NcN} = @SVector [C.constr1.link.data.id]
@inline linkids(C::Combined{T,Nc,N,Nc²,NcN,2}) where {T,Nc,N,Nc²,NcN} = @SVector [C.constr1.link1.data.id, C.constr1.link2.data.id]
