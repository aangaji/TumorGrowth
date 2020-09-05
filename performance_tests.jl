using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise, BenchmarkTools

@time using TumorGrowth

#########################
### Performance tests ###
#########################

using StaticArrays

@benchmark pos = [ rand(3) for _=1:10000]
@benchmark posstatic = [@SVector rand(3) for _=1:10000]


pos = [rand(3) for _=1:10000]
posstatic = [@SVector rand(3) for _=1:10000]

pos = [(p .- 0.5)*50. for p in pos]
posstatic = [(p .- 0.5)*50. for p in posstatic]


N = length(pos)

@benchmark TumorGrowth.pos2box.([1.,1.,2.])

@benchmark [TumorGrowth.pos2box.(p) for p in pos]
@benchmark [TumorGrowth.pos2box.(p) for p in posstatic]

cellbox = [TumorGrowth.pos2box.(p) for p in posstatic]

pos2box(x::Float64) = Int.(ceil.( ( x-1, x+1 )./4))
pos2box_2(x::Float64) = (Int(ceil((x-1)/4)), Int(ceil((x+1)/4)))

@benchmark [pos2box.(p) for p in posstatic]
@benchmark [pos2box_2.(p) for p in posstatic]

root = 4

TumorGrowth.find_neighbors(cellbox, root)

@benchmark (1:N)[TumorGrowth.find_neighbors(cellbox, root)]

@benchmark TumorGrowth.find_neighbors(cellbox, root) |> findall

iter = cellbox[root] .|> first |> enumerate

@benchmark (n!=root && all(any(ir-2 .<i[k].<ir+2) for (k,ir) in iter) for (n,i) in enumerate(cellbox)) |> collect

function find_neighbors(cellbox, root)
    iter = cellbox[root] .|> first |> enumerate
    return collect( n!=root && all(any(ir-2 .<i[k].<ir+2) for (k,ir) in iter) for (n,i) in enumerate(cellbox))
end

function find_neighbors_new(cellbox, root)
    iter = cellbox[root] .|> first |> enumerate
	bools = isneighbor.(cellbox; rootiter=iter)
	bools[root] = false
	return findall(bools)
end
isneighbor(singlecellbox; rootiter) = all(any(ir-2 .<singlecellbox[k].<ir+2) for (k,ir) in rootiter)

iter = cellbox[root] .|> first |> enumerate

@benchmark isneighbor.(cellbox; rootiter=iter)

@benchmark find_neighbors(cellbox, root)

@benchmark find_neighbors_new(cellbox, root)

@benchmark birth_death_pushing(300; b=1.,d=0.,mu=0.,dim=3)

@benchmark birth_death_pushing(500; b=1.,d=0.,mu=0.,dim=3)

birth_death_pushing(20000; b=1.,d=0.,mu=0.,dim=3)

@code_warntype birth_death_pushing(300; b=1.,d=0.,mu=0.,dim=3)

birth_death_pushing(2000; b=1.,d=0.,mu=0.,dim=3)[2]

typeof(MVector(1.,1.,1.))
MVector{d, Float64} where d

@MVector zeros(Float64, 3)

nothing


# struct Neighbors
# 	root::Int
# 	root_itr
# 	cellbox
# end
#
# function Base.iterate(it::Neighbors, n=1)
# 	if n > length(it)
# 		return nothing
# 	else
# 		return n!=it.root && all(any(ir-2 .<it.cellbox[n][k].<ir+2) for (k,ir) in iter), n+1
# 	end
# end
#
# Base.length(iter::Neighbors) = length(iter.cellbox)
#
# Base.eltype(iter::Neighbors) = Bool
#
# collect(Neighbors(root, iter, cellbox))
# for n in Neighbors(root, iter, cellbox)
# 	println(n)
# end
