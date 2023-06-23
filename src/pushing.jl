pos2box(x::Float64; w=1) = Int(fld(x, w)) + 500
pos2box(p::Vector{Float64}, dimv::Val; w=1) = CartesianIndex(ntuple(i -> pos2box(p[i]; w=w), dimv))

function find_neighbors(box2cell::Array{Int,dim}, cell::Cell, 
        exclude_cell::Val{true}; s=2, dimv=Val(dim)) where {dim}
    nh = pos2box(cell.position, dimv) .+ CartesianIndices(ntuple(i -> -s:s, dimv))
    return filter(i -> !iszero(i) && i != cell.index, box2cell[nh])
end
function find_neighbors(box2cell::Array{Int,dim}, cell::Cell,
        exclude_cell::Val{false}; s=2, dimv=Val(dim)) where {dim}
    nh = pos2box(cell.position, dimv) .+ CartesianIndices(ntuple(i -> -s:s, dimv))
    return filter(!iszero, box2cell[nh])
end
function find_neighbors(box2cell::Array{Int,dim}, cell::Cell,
        ; s=2, dimv=Val(dim)) where {dim}
    find_neighbors(box2cell, cell, Val(true); s=s, dimv)
end

function pushing!(tumor, root, box2cell, dimv::Val)
    queue = [root]
    while !isempty(queue)
        root = popfirst!(queue)
        c1 = tumor[root]
        r1 = c1.position

        for n in shuffle!(find_neighbors(box2cell, c1,Val(false); s=2))
            r2 = tumor[n].position
            d = r2.-r1
            if norm(d) < 2.
                
                box2cell[pos2box(r2, dimv)] = 0
                r2 .+= d*(2.2/norm(d) - 1)
                
                push!(queue, n)
            end
        end
        box2cell[pos2box(r1, dimv)] = root
    end
end

# function pushing_recursive!(tumor, root, cellbox, dimv::Val)
#     r1 = tumor[root].position

#     N=length(tumor)
#     for n in shuffle!(findall(find_neighbors(cellbox, root)))
#         r2 = tumor[n].position
#         d = r2.-r1
#         if norm(d) < 2.
#             r2 .+= d*(2.2/norm(d) - 1)
#             cellbox[n] = pos2box(r2, dimv)

#             pushing_recursive!(tumor, n, cellbox)
#         end
#     end
# end
