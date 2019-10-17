using VanDeVusseUnderUncertainty, PolyChaos, QuadGK

# auxiliary functions
function createTuples(degree::Int)
	tuples = []
	for i in 0:degree
		for j in 0:i
			push!(tuples,(i,j))
		end
	end
	tuples
end

function createCentralIndices(degree::Int)
	central = [1]
	for i in 2:degree+1
		push!(central, central[end]+i)
	end
	central
end

# main
degree = 4
α, β, w = [2, 4], [4.5, 1.5], [0.3, 0.7]

op_k = constructBasis(degree,α,β,w)

showbasis(op_k)

# the following line takes some time
# especially for cases that should be zero
values_of_integrals = [ quadgk(t -> evaluate(i,t,op_k)*evaluate(j,t,op_k)*op_k.measure.w(t), op_k.measure.dom...)[1] for (i,j) in createTuples(degree) ]
# extract values that should be non-zero
@show non_zero_entries = [ values_of_integrals[i] for i in createCentralIndices(degree) ]
# extract values that should be zero
@show zero_entries = [ values_of_integrals[i] for i in 1:length(createTuples(degree)) if i ∉ createCentralIndices(degree) ]


