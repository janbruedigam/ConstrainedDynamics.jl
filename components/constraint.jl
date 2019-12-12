abstract type Constraint{T,Nc,Nc²,Nl} <: Node{T,Nc} end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, C::Constraint)
    summary(io, C); println(io)
    print(io, "\nConnected links: ")
    show(io, mime, linkids(C))
end

Base.show(io::IO, C::Constraint) = summary(io, C)
