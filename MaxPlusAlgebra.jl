module MaxPlusAlgebra

export ℝᵋ, ε, ⊕, ⊙

struct ℝᵋ{T <: Real} <: Real
    λ::T
end

ℝᵋ(x::ℝᵋ) = ℝᵋ(x.λ)

const global ε = ℝᵋ(typemin(Float64))

Base.convert(::Type{Real}, x::ℝᵋ) = x.λ
Base.promote_rule(::Type{ℝᵋ{T}}, x::Type{T}) where {T<:Number} = ℝᵋ{T}(x)
Base.typemin(::Type{ℝᵋ{T}}) where {T<:Number} = ℝᵋ(typemin(T))

⊙(x::ℝᵋ, y::ℝᵋ) = ((x == ε || y == ε) && return ε; ℝᵋ(x.λ+y.λ))
⊙(x::Real, y::Real) = ℝᵋ(x+y)
⊙(x::Number, y::Number) = ⊙(promote(x,y)...)


⊕(x::ℝᵋ, y::ℝᵋ) = ℝᵋ(max(x.λ, y.λ))
⊕(x::Real, y::Real) = ℝᵋ(max(x,y))
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
