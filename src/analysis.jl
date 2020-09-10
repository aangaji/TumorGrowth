export mutation_freqs, clone, clones

function mutation_freqs(tumor)
    allmuts = vcat(tumor.mutations...)
    bins = sort(unique(allmuts)) |> u -> vcat(u, u[end]+1)
    hist = fit(Histogram, allmuts, bins, closed=:left)
    return hist.weights ./ size(tumor, 1)
end


function find_mut(tumor, mut)
    isnothing(mut) && return 1
    for muts in tumor.mutations
        i = findfirst(isequal(mut), muts)
        isnothing(i) || return i
    end
end

# clone returns all cells with "mut" as their "sub"th mutation. Cells with less than "sub" mutations are returned for mut=nothing

clone(tumor, mut; sub = find_mut(tumor, mut) ) = filter(r -> isnothing(sub) || length(r.mutations)<sub ? isnothing(mut) : getindex(r.mutations, sub) == mut, tumor )

# clones returns all subclones with unique mutations as "sub"th mutation and one subclone with less than "sub" mutations (clonal or no mutations)

function clones(tumor; sub=1)

    out = Vector{DataFrame}()
    !isnothing(findfirst(length.(tumor.mutations).<sub)) && push!(out, clone(tumor, nothing, sub=sub))

    muts = getindex.(filter(r -> length(r.mutations)>=sub, tumor).mutations, sub) |> unique!

    for mut in muts
        push!(out, clone(tumor, mut, sub=sub))
    end

    length(out)==1 && return clones(tumor; sub=sub+1)
    return out
end
