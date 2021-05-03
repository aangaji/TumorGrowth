export mutation_freqs, stochastic_sequencing, multi_region_sequencing, reduced_μ, reduced_μ!, clone, clones, haplotypes, clones_by_mutations

function mutation_freqs(tumor)
    allmuts = vcat(tumor.mutations...)
    isempty(allmuts) && return DataFrame(mutation = Int[],
                                        reads = Int[],
                                        coverage = Int[],
                                        frequency = Float64[])
    muts = union(allmuts)
    bins = vcat(muts, muts[end]+1)
    hist = fit(Histogram, allmuts, bins, closed=:left)

    DataFrame(mutation = muts,
            reads = hist.weights,
            coverage = fill(nrow(tumor), length(muts)),
            frequency = hist.weights ./ nrow(tumor)
            )
end

function stochastic_sequencing(tumor; readdepth = 100)
    seq_results = mutation_freqs(tumor)
    seq_results.reads = [ rand( Poisson(readdepth*f) ) for f = seq_results.frequency ]
    seq_results.coverage = fill(readdepth, nrow(seq_results))
    seq_results.frequency = seq_results.reads ./ seq_results.coverage
    seq_results
end


function multi_region_sequencing(tumor; n=0, a=0., cells_per_sample=0, sample_r=a/2, res=0.0, stochastic = false, readdepth = 100)
    lattice, samples = multi_region_sampling(tumor; n=n, a=a, cells_per_sample=cells_per_sample, sample_r=sample_r)
    seq_results = filter!.(c->c.frequency > res, stochastic ? stochastic_sequencing.(samples; readdepth=readdepth) : mutation_freqs.(samples))
    sampletumor = DataFrame(
        index = 1:length(samples),
        n = nrow.(samples),
        position = lattice,
        mutations = getproperty.(seq_results, :mutation),
        frequencies = getproperty.(seq_results, :frequency),
        reads = getproperty.(seq_results, :reads),
        coverages = getproperty.(seq_results, :coverage)
        )

    return samples, sampletumor
end


function reduced_μ!(tumor, x)
    tumor |> mutation_freqs |>
        seq -> filter!(_-> rand()<=x, seq) |>
        seq_red -> begin
            for muts in tumor.mutations
                filter!(m-> m in seq_red.mutation, muts)
            end
        end
    tumor
end
reduced_μ(tumor, x) = reduced_μ!(deepcopy(tumor), x)

function haplotypes(tumor; res=0.0)
    type_muts = unique(tumor.mutations)
    filter!(!isempty, type_muts)
    types = type_muts .|> m-> filter(c-> c.mutations == m , tumor)
    select = nrow.(types)./nrow(tumor) .>= res
    return types[select], type_muts[select]
end

function clones_by_mutations(tumor; res=0.0)
    isempty(tumor) && return (Vector{typeof(tumor)}(), tumor.mutations)
    mutations = union(tumor.mutations...)

    types = mutations .|> m -> filter(r -> m in r.mutations, tumor)
    select = nrow.(types)./nrow(tumor) .>= res
    return types[select], mutations[select]
end

############################
# The following methods are intended for plotting of tumors with mutations
# sorted in order of occurence

function find_mut(tumor, mut)
    isnothing(mut) && return 1
    for muts in tumor.mutations
        i = findfirst(isequal(mut), muts)
        isnothing(i) || return i
    end
end

# clone returns all cells with "mut" as their "sub"th mutation. Cells with less than "sub" mutations are returned for mut=nothing

function clone(tumor, mut; sub = find_mut(tumor, mut), show_warning = true )
    show_warning && @warn "This method assumes mutations to be ordered by occurence!"
    filter(r -> isnothing(sub) || length(r.mutations)<sub ? isnothing(mut) : getindex(r.mutations, sub) == mut, tumor )
end


# clones returns all subclones with unique mutations as "sub"th mutation and one subclone with less than "sub" mutations (clonal or no mutations)

function clones(tumor; sub=1, autodepth=true, show_warning = true)

    out = Vector{DataFrame}()
    !isnothing(findfirst(l-> l<sub, length.(tumor.mutations))) && push!(out, clone(tumor, nothing, sub=sub, show_warning=false))

    muts = getindex.(filter(r -> length(r.mutations)>=sub, tumor).mutations, sub) |> unique!

    autodepth && (isempty(out) && length(muts)<2) && return clones(tumor; sub=sub+1)
    show_warning && @warn "This method assumes mutations to be ordered by occurence!"

    for mut in muts
        push!(out, clone(tumor, mut; sub=sub, show_warning = false))
    end
    return out
end
