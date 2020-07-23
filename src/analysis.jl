export allele_population

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


### VAF distribution
# nBins = 50
# hist = fit(Histogram, f_select, nbins=nBins, closed=:left)
# bar(midpoints(hist.edges[1]), hist.weights, xlims=(0.,1.1), xlabel=:f)
