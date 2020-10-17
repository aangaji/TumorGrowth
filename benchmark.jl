using Pkg; Pkg.activate(pwd()); Pkg.instantiate()

using BenchmarkTools

@time using TumorGrowth
    # 83.128724 seconds (40.55 M allocations: 2.109 GiB, 1.01% gc time)
    # 27.123793 seconds (58.74 M allocations: 3.013 GiB, 4.57% gc time)
    # 50.520389 seconds (51.40 M allocations: 2.660 GiB, 2.21% gc time)

@benchmark birth_death_pushing(2000; b=0.69, d=0.1, mu=0.3, ρc=1., dim=2, seed = 1010) seconds=120

# BenchmarkTools.Trial:
#   memory estimate:  3.40 GiB
#   allocs estimate:  141194144
#   --------------
#   minimum time:     5.266 s (3.87% GC)
#   median time:      5.344 s (4.55% GC)
#   mean time:        5.333 s (4.43% GC)
#   maximum time:     5.381 s (4.73% GC)
#   --------------
#   samples:          23
#   evals/sample:     1

@benchmark birth_death_pushing(2000; b=0.69, d=0.1, mu=0.3, ρc=1., dim=2, seed = 8476) seconds=120

# BenchmarkTools.Trial:
#   memory estimate:  3.35 GiB
#   allocs estimate:  139509671
#   --------------
#   minimum time:     5.209 s (3.56% GC)
#   median time:      5.265 s (4.48% GC)
#   mean time:        5.269 s (4.34% GC)
#   maximum time:     5.330 s (4.07% GC)
#   --------------
#   samples:          23
#   evals/sample:     1

@benchmark birth_death_pushing(2000; b=0.69, d=0.1, mu=0.3, ρc=20., dim=2, seed = 1010) seconds=120

# BenchmarkTools.Trial:
#   memory estimate:  1.81 GiB
#   allocs estimate:  8996400
#   --------------
#   minimum time:     2.453 s (6.21% GC)
#   median time:      2.492 s (6.46% GC)
#   mean time:        2.491 s (6.41% GC)
#   maximum time:     2.528 s (6.58% GC)
#   --------------
#   samples:          49
#   evals/sample:     1

@benchmark birth_death_pushing(2000; b=0.69, d=0.1, mu=0.3, ρc=20., dim=2, seed = 8476) seconds=120

# BenchmarkTools.Trial:
#   memory estimate:  1.89 GiB
#   allocs estimate:  9252040
#   --------------
#   minimum time:     2.608 s (5.79% GC)
#   median time:      2.658 s (6.31% GC)
#   mean time:        2.656 s (6.25% GC)
#   maximum time:     2.691 s (6.67% GC)
#   --------------
#   samples:          46
#   evals/sample:     1

@benchmark birth_death_pushing(5000; b=0.69, d=0.1, mu=0.3, ρc=7., dim=3, seed = 1010) seconds=120

# BenchmarkTools.Trial:
#   memory estimate:  4.26 GiB
#   allocs estimate:  160396621
#   --------------
#   minimum time:     7.610 s (4.16% GC)
#   median time:      7.672 s (5.22% GC)
#   mean time:        7.686 s (5.07% GC)
#   maximum time:     7.758 s (5.59% GC)
#   --------------
#   samples:          16
#   evals/sample:     1

@benchmark birth_death_pushing(5000; b=0.69, d=0.1, mu=0.3, ρc=7., dim=3, seed = 8476) seconds=120

# BenchmarkTools.Trial:
#   memory estimate:  4.35 GiB
#   allocs estimate:  162446164
#   --------------
#   minimum time:     7.848 s (4.21% GC)
#   median time:      7.944 s (4.81% GC)
#   mean time:        7.939 s (4.74% GC)
#   maximum time:     8.045 s (5.26% GC)
#   --------------
#   samples:          16
#   evals/sample:     1

@benchmark birth_death_pushing(5000; b=0.69, d=0.1, mu=0.3, ρc=30., dim=3, seed = 1010) seconds=120

# BenchmarkTools.Trial:
#   memory estimate:  2.90 GiB
#   allocs estimate:  51563497
#   --------------
#   minimum time:     7.282 s (3.60% GC)
#   median time:      7.359 s (4.08% GC)
#   mean time:        7.351 s (4.00% GC)
#   maximum time:     7.414 s (4.14% GC)
#   --------------
#   samples:          17
#   evals/sample:     1

@benchmark birth_death_pushing(5000; b=0.69, d=0.1, mu=0.3, ρc=30., dim=3, seed = 8476) seconds=120

# BenchmarkTools.Trial:
#   memory estimate:  2.92 GiB
#   allocs estimate:  51053037
#   --------------
#   minimum time:     7.475 s (3.52% GC)
#   median time:      7.523 s (3.61% GC)
#   mean time:        7.522 s (3.58% GC)
#   maximum time:     7.578 s (3.78% GC)
#   --------------
#   samples:          16
#   evals/sample:     1

using SharedArrays, Distributed
addprocs(5)
workers()
@everywhere println("hello")
#rmprocs(workers())


@everywhere using Pkg
@everywhere Pkg.activate(pwd())
@everywhere using BenchmarkTools
@everywhere using TumorGrowth

using Plots
fig = plot()
for d = [0.0, 0.1, 0.2]
    ρ = 6.
    Nrange = 1000:1000:10000
    times = SharedVector(zeros(Float64,length(Nrange)))
    for (i,N) in enumerate(Nrange)
        println("d $d, N $N")
        times[i] = @distributed (+) for _=1:5
            repeat = true
            t = 0.
            while repeat
                try
                    t = @elapsed birth_death_pushing(N; b=0.69, d=d, mu=0.3, ρc=ρ, dim=2, seed = rand(1000:9999))
                    repeat = false
                catch e
                end
            end
            t/5
        end
    end
    plot!(fig, Nrange, times, lab="d=$d, ρc=$ρ", marker=:o)
end
plot!(xlabel=:N, ylabel="time/s", legend=:topleft)
#savefig(fig, "runtimes3d.png")
