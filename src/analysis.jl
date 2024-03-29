export mutation_freqs, sampletumor_mfreqs, stochastic_sequencing, multi_region_sequencing, reduced_mu, reduced_mu!, clone, clones, haplotypes, clones_by_mutations

function mutation_freqs(tumor)
    all(isempty.(tumor.mutations)) && return DataFrame(mutation = Int[],
                                    reads = Int[],
                                    coverage = Int[],
                                    frequency = Float64[])
    vcat(tumor.mutations...) |> countmap |> collect |> sort |> list->
            DataFrame(
                    mutation=getindex.(list,1),
                    reads = getindex.(list, 2),
                    coverage = fill(nrow(tumor), length(list)),
                    frequency = getindex.(list, 2) ./ nrow(tumor)
                    )
end

function sampletumor_mfreqs(sampletumor)
    mutations = vcat(sampletumor.mutations...) |> unique! |> sort!
    m_ind = Dict(mutations .=> 1:length(mutations))
    freqs = zeros(length(mutations))
    for row in eachrow(sampletumor)
        for (m,f) in zip(row.mutations, row.frequencies)
            freqs[m_ind[m]] += f
        end
    end
    freqs./=nrow(sampletumor)
    return DataFrame(mutation = mutations, frequency=freqs)
end

function stochastic_sequencing(tumor; readdepth)
    seq_results = mutation_freqs(tumor)
    seq_results.reads = [ rand( Binomial(readdepth,f) ) for f = seq_results.frequency ]
    seq_results.coverage = fill(readdepth, nrow(seq_results))
    seq_results.frequency = seq_results.reads ./ seq_results.coverage
    seq_results
end


function multi_region_sequencing(tumor; n=0, a=0., cells_per_sample=0, sample_r=a/2, res=0.0, stochastic = false, readdepth = 0)
    lattice, samples, sample_r = multi_region_sampling(tumor; n=n, a=a, cells_per_sample=cells_per_sample, sample_r=sample_r)
    seq_results = filter!.(c->c.frequency > res, stochastic ? stochastic_sequencing.(samples; readdepth=readdepth) : mutation_freqs.(samples))
    sampletumor = DataFrame(
        index = 1:length(samples),
        n = nrow.(samples),
        sample_r = fill(sample_r, length(samples)),
        position = lattice,
        mutations = getproperty.(seq_results, :mutation),
        frequencies = getproperty.(seq_results, :frequency),
        reads = getproperty.(seq_results, :reads),
        coverages = getproperty.(seq_results, :coverage)
        )

    return (samples=samples, sampletumor=sampletumor)
end


function reduced_mu!(tumor, x)
    reduced = vcat(tumor.mutations...) |> unique! |> sort! |>
        muts -> filter!(muts) do m
            rand()<=x
        end
    filter!.(in(reduced), tumor.mutations)
    return tumor
end

reduced_mu(tumor, x) = reduced_mu!(deepcopy(tumor), x)

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
