export mutation_freqs, multi_region_sequencing, reduced_μ, reduced_μ!, clone, clones, haplotypes, clones_by_mutations

function mutation_freqs(tumor)
    allmuts = vcat(tumor.mutations...)
    isempty(allmuts) && return DataFrame(mutation = Int[], frequency = Float64[])
    muts = allmuts |> unique |> sort
    bins = vcat(muts, muts[end]+1)
    hist = fit(Histogram, allmuts, bins, closed=:left)
    return DataFrame(mutation = muts , frequency = hist.weights ./ size(tumor, 1))
end

function multi_region_sequencing(tumor; n=0, a=0., cells_per_sample=0, sample_r=a/2, res=0.0)
    lattice, samples = multi_region_sampling(tumor; n=n, a=a, cells_per_sample=cells_per_sample, sample_r=sample_r)
    seq_results = filter!.(c->c.frequency > res, mutation_freqs.(samples))
    sampletumor = DataFrame(index = 1:length(samples), position = lattice, mutations = getproperty.(seq_results, :mutation), frequencies = getproperty.(seq_results, :frequency))

    return samples, sampletumor
end


function reduced_μ!(tumor, x)
    tumor |> mutation_freqs |>
        seq -> filter!(_-> rand()<x, seq) |>
        seq_red -> begin
            for muts in tumor.mutations
                filter!(m-> m in seq_red.mutation, muts)
            end
        end
    tumor
end
reduced_μ(tumor, x) = reduced_μ!(deepcopy(tumor), x)


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

function haplotypes(tumor; res=0.0)
    mutations = filter!(c-> c.frequency > res, mutation_freqs(tumor)).mutation

    types = mutations .|> m -> filter(r -> !isempty(r.mutations) && last(r.mutations) == m, tumor)
    filter!(t -> !isempty(t), types), getproperty.(first.(types),:mutations)
end

function clones_by_mutations(tumor; res=0.0)
    mutations = filter!(c-> c.frequency > res, mutation_freqs(tumor)).mutation

    types = mutations .|> m -> filter(r -> m in r.mutations, tumor)
    types, mutations
end
