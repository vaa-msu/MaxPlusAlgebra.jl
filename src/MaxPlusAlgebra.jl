module MaxPlusAlgebra

export ℝᵋ, ε, ⊕, ⊙

struct ℝᵋ <: Real
    value::Float64
end

ℝᵋ(x::ℝᵋ) = ℝᵋ(x.value)

Base.max(x::ℝᵋ, y::ℝᵋ) = ℝᵋ(max(x.value, y.value))
Base.:+(x::ℝᵋ, y::ℝᵋ) = ℝᵋ(x.value + y.value)
Base.typemin(::Type{ℝᵋ}) = ℝᵋ(typemin(Float64))

const global ε = typemin(ℝᵋ)

Base.promote_rule(::Type{ℝᵋ}, ::Type{T}) where {T<:Number} = ℝᵋ

⊙(x::T, y::T) where {T<:Real} = x + y
⊙(x::Number, y::Number) = ⊙(promote(x,y)...)

⊕(x::T, y::T) where {T<:Real} = max(x,y)
⊕(x::Number, y::Number) = ⊕(promote(x,y)...)


ℝᵋ(A::Array) = map(ℝᵋ, A)

function ⊡(A::Matrix{T}, B::Matrix{T}, ⊞ = +, ⊠ = *, ∅ = 0) where T<:Number
    @assert size(A)[end] == size(B)[1]
    S = Matrix{T}(undef, size(A)[1], size(B)[end])
    for i = 1:size(A)[1], j = 1:size(B)[end]
        S[i,j] = ∅
        for k = 1:size(A)[end]
            S[i,j] = ⊞(S[i,j], ⊠(A[i,k], B[k,j]))
        end
    end
    return S
end

⊙(A::Matrix{T}, B::Matrix{T}) where T<:Number = ⊡(A, B, ⊕, ⊙, typemin(T))

end
