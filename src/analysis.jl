export allele_population, mutation_freqs, clone, clones

function allele_population(mutations; skip=Int64[])
    muts = Array{Int64}(undef,0)
    for c in mutations
        for m in c
            m in skip || push!(muts, m)
        end
    end
    return muts
end

function mutation_freqs(mutations; skip=Int[])
    muts = allele_population(mutations; skip=skip)
    mut_counts = [count(isequal(m), muts) for m in sort(unique(muts)) ]
    return mut_counts./ length(mutations)
end



function find_mut(tumor, mut)
    isnothing(mut) && return 1
    for muts in tumor.mutations
        i = findfirst(isequal(mut), muts)
        isnothing(i) || return i
    end
end

# clone returns all cells with "mut" as their "sub"th mutation. Cells with less than "sub" mutations are returned for mut=nothing

clone(tumor, mut; sub = find_mut(tumor, mut) ) = filter(r -> length(r.mutations)<sub ? isnothing(mut) : getindex(r.mutations, sub) == mut, tumor )

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
