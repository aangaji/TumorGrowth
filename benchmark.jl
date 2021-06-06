using Pkg; Pkg.activate(pwd()); Pkg.instantiate()

using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 300

@time using TumorGrowth
    # 83.128724 seconds (40.55 M allocations: 2.109 GiB, 1.01% gc time)
    # 27.123793 seconds (58.74 M allocations: 3.013 GiB, 4.57% gc time)
    # 50.520389 seconds (51.40 M allocations: 2.660 GiB, 2.21% gc time)

function rerun(until; params...)
    try
        out = birth_death_pushing(until; params...)
        return out
    catch e
        e == ErrorException("Tumor died") || rethrow()
        return rerun(until; params... )
    end
end

@benchmark rerun(2000; b=0.69, d=0.1, μ=0.3, ρ=1., dim=2)

BenchmarkTools.Trial:
  memory estimate:  632.83 MiB
  allocs estimate:  4827004
  --------------
  minimum time:     1.212 s (10.21% GC)
  median time:      1.553 s (12.36% GC)
  mean time:        1.550 s (12.53% GC)
  maximum time:     2.759 s (41.91% GC)
  --------------
  samples:          194
  evals/sample:     1

@benchmark rerun(2000; b=0.69, d=0.1, μ=0.3, ρ=20., dim=2)

BenchmarkTools.Trial:
  memory estimate:  1.72 GiB
  allocs estimate:  4531734
  --------------
  minimum time:     2.046 s (8.06% GC)
  median time:      2.232 s (9.20% GC)
  mean time:        2.230 s (9.17% GC)
  maximum time:     2.654 s (9.13% GC)
  --------------
  samples:          135
  evals/sample:     1

@benchmark rerun(5000; b=0.69, d=0.1, μ=0.3, ρ=7., dim=3)

BenchmarkTools.Trial:
  memory estimate:  1.39 GiB
  allocs estimate:  16098472
  --------------
  minimum time:     3.047 s (5.51% GC)
  median time:      3.450 s (6.18% GC)
  mean time:        3.471 s (6.27% GC)
  maximum time:     4.075 s (16.21% GC)
  --------------
  samples:          87
  evals/sample:     1

@benchmark rerun(5000; b=0.69, d=0.1, μ=0.3, ρ=30., dim=3)

BenchmarkTools.Trial:
  memory estimate:  2.11 GiB
  allocs estimate:  13781166
  --------------
  minimum time:     5.545 s (4.07% GC)
  median time:      5.770 s (3.99% GC)
  mean time:        5.766 s (3.97% GC)
  maximum time:     5.999 s (4.06% GC)
  --------------
  samples:          53
  evals/sample:     1



using SharedArrays, Distributed
addprocs(5)
nworkers()
@everywhere println("hello")
# rmprocs(workers())

@everywhere begin
    using Pkg
    Pkg.activate(pwd())
    using BenchmarkTools
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 180
    using TumorGrowth
end

@everywhere function rerun(until; params...)
    try
        out = birth_death_pushing(until; params...)
        return out
    catch e
        e == ErrorException("Tumor died") || rethrow()
        return rerun(until; params... )
    end
end

using Plots

fig = plot(xlabel=:N, ylabel="time/s", legend=:topleft)
b, μ = 1., 0.3
dim, ρ = 3, 6.

for d = [0.0, 0.1, 0.2, 0.3]
    Nrange = 1000:1000:10000
    times = SharedVector(zeros(Float64,length(Nrange)))
    @sync @distributed for i=1:length(Nrange)
        times[i] = @belapsed rerun($(Nrange[i]); b=$b, d=$d, μ=$μ, ρ=$ρ, dim=$dim, showprogress=false)
        println(myid(),": d $d, N $(Nrange[i])")
    end
    plot!(fig, Nrange, times, lab="d=$d, ρ=$ρ", marker=:o)
    display(fig)
end
#savefig(fig, "examples/runtimes3d.png")
