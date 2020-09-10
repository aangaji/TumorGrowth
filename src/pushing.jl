pos2box(x::Float64) = Int.(ceil.( ( x-1, x+1 )./4))

# function find_neighbors(cellbox, root)
#     iter = cellbox[root] .|> first |> enumerate
#     return collect( n!=root && all(any(ir-2 .<i[k].<ir+2) for (k,ir) in iter) for (n,i) in enumerate(cellbox))
# end

function find_neighbors(cellbox, root)
    iter = cellbox[root] .|> first |> enumerate
	bools = isneighbor.(cellbox; rootiter=iter)
	bools[root] = false
	return bools
end
isneighbor(singlecellbox; rootiter) = all(any(ir-2 .<singlecellbox[k].<ir+2) for (k,ir) in rootiter)

function pushing!(tumor, root, cellbox)
    queue = [root]
    N = length(tumor)
    while !isempty(queue)
        root = popfirst!(queue)
        r1 = tumor[root].position

        for n in shuffle!(findall(find_neighbors(cellbox, root)))
            r2 = tumor[n].position
            d = r2.-r1
            if norm(d) < 2.
                r2 .+= d*(2.2/norm(d) - 1)#/2
                #r1 .-= d*(2.2/norm(d) - 1)/2
                cellbox[n] = pos2box.(r2)
                #cellbox[root] = pos2box.(r1)
                push!(queue, n)
            end
        end
    end
end

function pushing_recursive!(tumor, root, cellbox)
    r1 = tumor[root].position

    N=length(tumor)
    for n in shuffle!(findall(find_neighbors(cellbox, root)))
        r2 = tumor[n].position
        d = r2.-r1
        if norm(d) < 2.
            r2 .+= d*(2.2/norm(d) - 1)
            cellbox[n] = pos2box.(r2)

            pushing_recursive!(tumor, n, cellbox)
        end
    end
end
