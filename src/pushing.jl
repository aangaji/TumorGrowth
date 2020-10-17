pos2box(x::Float64; w=2) = floor(Int,x/w)
pos2box(p::Vector{Float64}, dimv::Val; w=2) = ntuple(i->pos2box(p[i];w=w), dimv)

function find_neighbors(cellbox, root::Int; s=2)
	bools = isneighbor.(cellbox; rootbox=cellbox[root], s=s)
	bools[root] = false
	return bools
end
isneighbor(neighbox::NTuple{dim,Int}; rootbox::NTuple{dim,Int}, s=2) where dim = all(-s < neighbox[i] - rootbox[i] < s for i=1:dim)

function pushing!(tumor, root, cellbox, dimv::Val)
    queue = [root]
    while !isempty(queue)
        root = popfirst!(queue)
        r1 = tumor[root].position

        for n in shuffle!(findall(find_neighbors(cellbox, root; s=2)))
            r2 = tumor[n].position
            d = r2.-r1
            if norm(d) < 2.
                r2 .+= d*(2.2/norm(d) - 1)#/2
                #r1 .-= d*(2.2/norm(d) - 1)/2
                cellbox[n] = pos2box(r2, dimv)
                #cellbox[root] = pos2box(r1;dimv=Val(dim))
                push!(queue, n)
				#push!(queue, n)
            end
        end
    end
end

function pushing_recursive!(tumor, root, cellbox, dimv::Val)
    r1 = tumor[root].position

    N=length(tumor)
    for n in shuffle!(findall(find_neighbors(cellbox, root)))
        r2 = tumor[n].position
        d = r2.-r1
        if norm(d) < 2.
            r2 .+= d*(2.2/norm(d) - 1)
            cellbox[n] = pos2box(r2, dimv)

            pushing_recursive!(tumor, n, cellbox)
        end
    end
end
