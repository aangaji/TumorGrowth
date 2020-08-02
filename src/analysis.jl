export allele_population, mutation_freqs

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
    mut_ids = sort(unique(muts))
    mut_counts = [count(isequal(m), muts) for m in mut_ids ]
    return mut_ids, mut_counts./ length(mutations)
end
