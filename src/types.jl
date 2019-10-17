export 	AbstractLTISystem,
		AbstractCost,
		AbstractUncertainty,
		LTISystem,
		Cost,
		Uncertainty

abstract type AbstractLTISystem end
abstract type AbstractCost end
abstract type AbstractUncertainty end

struct LTISystem <: AbstractLTISystem
	A::Array{<:Real,3}
	B::Array{<:Real,3}
	x0::Matrix
	nx::Int
	nu::Int
	L::Int
	δ::Real
	
	function LTISystem(A::Array{<:Real,3}, B::Array{<:Real,3}, x0::Matrix, δ::Real)
		nx0, lx0 = size(x0)
		nxA, nxA_, lA = size(A)
		nxB, nuB, lB = size(B)
		nxA != nxA_ && throw(DomainError(size(A),"system matrix must be square"))
		!(lA == lB == lx0) && throw(DomainError((lA,lB,lx0),"incosistent PCE dimensions"))
		!(nxA == nxB == nx0) && throw(DomainError((nxA,nxB,nx0),"incompatible sizes"))
		δ <= 0 && throw(DomainError(δ,"sampling time must be positive"))

		new(A, B, x0, nxA, nuB, lA, δ)
	end
end

struct Cost <: AbstractCost
	Q::Matrix
	R::Matrix

	function Cost(Q::Matrix, R::Matrix)
		!_issquare(Q) && throw(DomainError(size(Q),"cost matrix Q must be square"))
		!_issquare(R) && throw(DomainError(size(R),"cost matrix R must be square"))
		eigmin(Q) < 0 && throw(DomainError(eigmin(Q),"cost matrix Q must be positive semi-definite"))
		eigmin(R) <= 0 && throw(DomainError(eigmin(R),"cost matrix R must be positive definite"))

		new(Q,R)
	end
end

struct Uncertainty <: AbstractUncertainty
	op::AbstractOrthoPoly
	T2::Tensor
	T3::Tensor

	function  Uncertainty(mop::AbstractOrthoPoly)
		new(mop, Tensor(2, mop), Tensor(3, mop))
	end
end



