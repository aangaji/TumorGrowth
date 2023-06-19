pos2box(x::Float64; w=1) = Int(fld(x, w)) + 500
pos2box(p::Vector{Float64}, dimv::Val; w=1) = CartesianIndex(ntuple(i -> pos2box(p[i]; w=w), Val(2)))

function find_neighbors(box2cell::Array{Int,dim}, cell2box, root::Int; s=2) where {dim}
    nh = cell2box[root] .+ CartesianIndices(ntuple(i -> -s:s, Val(2)))
    return filter(i -> !iszero(i) && i != root, box2cell[nh])
end

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
