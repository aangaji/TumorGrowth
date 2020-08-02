pos2box(x::Float64) = Int.(ceil.( ( x-1, x+1 )./4))

function find_neighbors(cellbox, root)
    iter = cellbox[root] .|> first |> enumerate
    neighbors = Int[]
    for (n,i) in enumerate(cellbox)
        n!=root && all(any(ir-2 .<i[k].<ir+2) for (k,ir) in iter) && push!(neighbors, n)
    end
    return neighbors
end

function pushing!(tumor, root, cellbox)
    queue = [root]
    while !isempty(queue)
        root = popfirst!(queue)
        r1 = tumor.position[root]

        for n in shuffle!(find_neighbors(cellbox, root))
            r2 = tumor.position[n]
            d = norm(r1.-r2)
            if d < 2.
                r2 .+= (r2.-r1)*(2.1/d - 1)
                cellbox[n] = pos2box.(r2)
                tumor.position[n] = r2
                push!(queue, n)
            end
        end
    end
end

function pushing_recursive!(tumor, root, cellbox)
    r1 = tumor.position[root]

    for n in shuffle(find_neighbors(cellbox, root))
        r2 = tumor.position[n]
        d = norm(r1.-r2)
        if d < 2.
            r2 .+= (r2.-r1)*(2.1/d - 1)
            cellbox[n] = pos2box.(r2)
            tumor.position[n] = r2
            pushing_recursive!(tumor, n, cellbox)
        end
    end
end
