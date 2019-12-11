mutable struct FillInData{T,N1,N2,N1N2}
    pid::Int64
    cid::Int64
    JL::SMatrix{N2,N1,T,N1N2}
    JU::SMatrix{N1,N2,T,N1N2}

    function FillInData{T,N1,N2}(pid,cid) where {T,N1,N2}
        JL = @SMatrix zeros(T,N2,N1)
        JU = @SMatrix zeros(T,N1,N2)
        new{T,N1,N2,N1*N2}(pid,cid,JL,JU)
    end
end

struct FillIn{T,N1,N2,N1N2}
    data::FillInData{T,N1,N2,N1N2}

    function FillIn{T,N1,N2}(pid,cid) where {T,N1,N2}
        data = FillInData{T,N1,N2}(pid,cid)
        new{T,N1,N2,N1*N2}(data)
    end
end
