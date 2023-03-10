"""
This file demonstrates most functionality of the package.
Load Revise if you want to make changes to the source files.
"""

using Pkg;
Pkg.activate(pwd());
using Revise

@time using TumorGrowth

#########################
######## SANDBOX ########
#########################


#########################
####### SIMULATION ######
#########################

@time simoutput = birth_death_pushing(3000; b=1.0, d=0.0, μ=0.1, ρ=Inf, dim=2, seed=1010, showprogress=false)
tumor = simoutput[:tumor]
mutations = simoutput[:mutations]

@time simoutput = nonspatial(1000000; b=1.0, d=0.0, μ=0.3, showprogress=true)

bumor = deepcopy(tumor)

birth_death_pushing!(bumor, length(tumor) + 1; b=0.69, d=0.0, μ=0.3, dim=2,
        t=simoutput[:time], cur_id=simoutput[:index], cur_mutation=simoutput[:mutation])

birth_death_pushing(12.0; b=0.69, d=0.0, μ=0.3, dim=3, seed=1010)

tumordf = tumor |> DataFrame
mutationsdf = mutations |> DataFrame

#########################
##### TUMOUR SAVING #####
#########################

using DataFrames, CSV

tumor = birth_death_pushing(5000; b=1.0, d=0.0, μ=0.3, ρ=6.0, dim=3, seed=1234)[:tumor]
CSV.write("examples/test_5000_3d.csv", DataFrame(tumor), delim='\t')
b = data_import("examples/test_5000_3d.csv")

########################
####### PLOTTING #######
########################

using Plots: palette

b = data_import("examples/test_2000_2d.csv")
@time plotting(b)
@time plotting_colored_mutations(b, colorpalette=palette(:tab20), shading=false, inline=false, show_warning=true, path="examples/test_2000_2d.png")

b = data_import("examples/test_5000_3d.csv")
@time plotting(b)
@time plotting_colored_mutations(b, shading=true, path="examples/test_5000_3d.png", size=(500, 500))

# translating and layouts

using GLMakie

tumor1 = birth_death_pushing(10000; b=1.0, d=0.0, μ=0.2, ρ=7.0, dim=3)[:tumor] |> DataFrame
tumor2 = birth_death_pushing(10000; b=1.0, d=0.0, μ=0.2, ρ=Inf, dim=3)[:tumor] |> DataFrame
tumor3 = birth_death_pushing(3000; b=1.0, d=0.0, μ=0.2, ρ=1.4, dim=2)[:tumor] |> DataFrame

scene = Scene(resolution=(1000, 1000), show_axis=false)
t1 = plotting_colored_mutations!(scene, tumor1;
        shading=true, colorpalette=cgrad(:reds)[0.4:0.1:1.0])
t2 = plotting_colored_mutations!(scene, tumor2;
        shading=true, colorpalette=cgrad(:blues)[0.4:0.1:1.0])[end]
t3 = plotting_colored_mutations!(scene, tumor3;
        shading=false, colorpalette=cgrad(:greens)[0.2:0.1:1.0])[end]
translate!(t2, Vec3f0(100, 0, 0))
translate!(t3, Vec3f0(50, 100, 0))
scene |> display

using GLMakie.MakieLayout

scene, layout = layoutscene(show_axis=false)
layout[1, 1:2] = ax1 = LScene(scene, camera=cam2d!, resolution=(500, 500))
layout[2, 1] = ax2 = LScene(scene, camera=cam3d!, raw=false, show_axis=false, resolution=(500, 500))
layout[2, 2] = ax3 = LScene(scene, camera=cam3d!, raw=false, show_axis=false, resolution=(500, 500))

t1 = plotting_colored_mutations!(ax2, tumor1;
        shading=true, colorpalette=cgrad(:reds)[0.4:0.1:1.0])
t2 = plotting_colored_mutations!(ax3, tumor2;
        shading=true, colorpalette=cgrad(:blues)[0.4:0.1:1.0])
# b = deepcopy(tumor3); b.position ./= 1.8
t3 = plotting_colored_mutations!(ax1, tumor3; markersize=1.1,
        shading=false, colorpalette=cgrad(:greens)[0.2:0.1:1.0])

scene |> display

# plotting sampletumors

@time tumor = birth_death_pushing(3000; b=1.0, d=0.0, μ=0.5, ρ=1.6, dim=2, showprogress=false)[:tumor] |> DataFrame
samples, sampletumor = multi_region_sequencing(tumor, n=50)

plotting_colored_mutations(sampletumor; markersize=sampletumor.sample_r, show_axis=true)
plotting_sampletumor_pies(sampletumor, first=30, colorpalette=TumorGrowth.distinguishable_colors(30))
plotting_sampletumor_pies(sampletumor, mutations=[1, 10, 55])

########################
####### Sampling #######
########################

b = data_import("examples/test_2000_2d.csv")

fig = b |> plotting_colored_mutations
cross_section(b; x=10.0, width=3.0) |> data -> plotting!(fig, data; color=:black)
cross_section(b; y=25.0, width=3.0) |> data -> plotting!(fig, data; color=:black)
radial_sample(b; r=30.0, width=3.0) |> data -> plotting!(fig, data; color=:black, path="examples/test_2000_2d.png")

fig = b |> plotting_colored_mutations
bulk(b; pos=(25.0, 0.0), box=(20, 30)) |> data -> plotting!(fig, data; color=:black)
punch(b; pos=(-20, -20), r=10) |> data -> plotting!(fig, data; color=:black)


b = data_import("examples/test_5000_3d.csv")
scene = b |> plotting_colored_mutations
cross_section(b; x=10.0, width=3.0) |> data -> plotting!(scene, data; shading=true, color=:black, path="examples/test_5000_3d.png")
cross_section(b; y=5.0, width=3.0) |> data -> plotting!(scene, data; color=:black)

display(scene)

fig = cross_section(b; x=10.0, width=6.0) |> plotting_colored_mutations
radial_sample(b; r=20.0, width=3.0) |> data -> cross_section(data; x=10.0, width=6.0) |> data -> plotting!(fig, data; color=:black)

########################
####### Analysis #######
########################

using Plots, StatsBase, Statistics

params = (b=1.0, d=0.2, μ=0.3, ρ=Inf, dim=3)
tumor = birth_death_pushing(4000; params...)[:tumor] |> DataFrame

using Plots, StatsBase, Statistics

function M!(fig, f; res=0.0, nBins=100, normalize=false, plotargs...)
        hist = fit(Histogram, 1 ./ f[res.<f], nbins=nBins, closed=:left)
        Mcounts = [sum(hist.weights[1:n]) for n = 1:length(hist.weights)]
        Mcounts = normalize ? Mcounts ./ last(Mcounts) : Mcounts
        Plots.plot!(fig, midpoints(hist.edges[1]), Mcounts; plotargs...)
end
M(f; res=0.0, nBins=100, normalize=false, plotargs...) = M!(Plots.plot(), f; res=res, nBins=nBins, normalize=normalize, plotargs...)

fig = M(mutation_freqs(tumor).frequency; res=0.0005, nBins=200,
        xlab="1/f", ylab="M(f)", lab="", marker=:o, ms=1.5)
β = params.μ * params.b / (params.b - params.d)
Plots.plot!(1:1000, n -> β * n, lab="")

f = filter(m -> m.reads > 1, stochastic_sequencing(tumor, readdepth=1000)).frequency
M!(fig, stochastic_sequencing(tumor, readdepth=2000).frequency; res=0.001, lab="", marker=:d, ms=1.5)

########################
######## clones ########
########################

tumor = data_import("examples/test_2000_2d.csv")

fig = plotting(tumor; color=:lightblue, shading=false)
clone(tumor, nothing) |> t -> plotting!(fig, t, color=:red)
clone(tumor, 2) |> t -> plotting!(fig, t, color=:blue)

fig = plotting_colored_mutations(tumor; colorpalette=palette(:tab20), shading=false)
sort(by=size, rev=true, clones(tumor))[1] |> t -> plotting!(fig, t, color=:black)

fig = plotting(tumor; color=:lightgrey, shading=false)
sort(by=size, rev=true, clones(tumor))[1] |> clones .|> t -> plotting_colored_mutations!(fig, t; colorpalette=palette(:tab10), shading=false, autodepth=true)
fig

#########################
###### time series ######
#########################

using DataFrames, GLMakie
using Plots: palette, distinguishable_colors


N, simparams... = (N=5000, b=1.0, d=0.0, μ=0.3, ρ=6.0, dim=3, seed=1234)
tumor = birth_death_pushing(N; simparams...)[:tumor] |> DataFrame

time_series = tumor_stepper(range(0.0, last(tumor.t_birth), length=50); simparams...)

last(time_series) == tumor

record_growth(time_series; path="temp.gif",
        frames=1, shading=true, show_axis=false,
        colorpalette=palette(:tab20)
)


#########################
###### density ~ b ######
#########################

using Plots

TumorGrowth.b_linear()
TumorGrowth.b_hill(6)
TumorGrowth.b_curve(1.0; bup=1.0, ρc=1.1)

d, ρ = 0.2, 1.5
@time tumor = birth_death_pushing(5000; b=0.69, d=d, μ=0.3, ρ=ρ, dim=2)[:tumor]
tumordf = tumor |> DataFrame

# power law / surface growth
Plots.plot(log.(tumordf.t_birth), log.(1:length(tumor)), legend=:none)
Plots.plot!(log.(tumordf.t_birth), t -> 2 * t, xlims=(0.0, 4.0))

# exponential / bulk
Plots.plot(tumordf.t_birth, log.(1:5000), legend=:none)
Plots.plot!(tumordf.t_birth, t -> (log(2) - d) * t)

Plots.plot(-7.0:0.01:7.0, r -> TumorGrowth.w(r; σ=3.0), fill=true)
Plots.plot!(-9.0:0.01:9.0, r -> TumorGrowth.w(r; σ=3.0), xticks=-7:2:7, legend=:none)

TumorGrowth.b_linear()
TumorGrowth.b_hill(6)
Plots.plot(0.0:0.01:9.0, ρ -> TumorGrowth.b_curve(ρ; bup=1.0, ρc=1.0), legend=:none)
Plots.vline!([1.0])

plotting_colored_mutations(tumordf, colorpalette=palette(:tab20), shading=false, inline=false)

Plots.histogram(tumordf.b, alpha=0.3, lw=0.2, lab="b")
Plots.vline!([median(tumordf.b)], lab="median", c=:blue)
Plots.vline!([d], lab="d", c=:black)

plotting(tumordf; color=tumordf.b ./ maximum(tumordf.b), colormap=cgrad(:reds), path="examples/surface_b.png")


#############################
### mulit region sampling ###
#############################

using Makie: meshscatter, meshscatter!, save, Point
using Statistics, DataFrames, LinearAlgebra, StatsBase

@time tumor = birth_death_pushing(5000; b=0.69, d=0.1, μ=0.3, ρ=1.8, dim=2, seed=1010)[:tumor] |> DataFrame

lattice, samples = multi_region_sampling(tumor; n=24, cells_per_sample=20);

length(lattice)

samples, sampletumor = multi_region_sequencing(tumor; n=100, cells_per_sample=20, stochastic=true)
length(samples)

mean(size.(samples, 1))



scene = plotting_colored_mutations(tumor)
for sample in samples
        plotting!(scene, sample; color=:black)
end
display(scene))

sample_r = 2 * sqrt(20 / π)
r = norm(std(tumor.position)) / 2
cm = mean(tumor.position)
density = size(punch(tumor; pos=cm, r=r), 1) / (π * r^2)
R = sqrt(size(tumor, 1) / (π * density))
# scene = plotting_colored_mutations(tumor)
a = minimum(pairwise(norm∘-,lattice)+I*100)
meshscatter!(scene, [Point{2}.(lattice)...], markersize=sample_r, scale_plot=false)
meshscatter!(scene, 0:0.2:360 .|> ϕ -> Point{2}((R-a/2) .* (cosd(ϕ), sind(ϕ)) .+ cm), markersize=1.0)
display(scene)

# save("multi_region_sampling.png", scene)

density = size(punch(tumor; pos=mean(tumor.position), r=20.0), 1) / (π * 20^2)
scene = plotting_colored_mutations(sampletumor; markersize=sqrt(20 / density / π), shading=true)
plotting_colored_mutations!(scene, tumor)

# save("multi_region_sequencing.png", scene)

##################################
## Reduce mutation artificially ##
##################################

b, d, μ = 0.69, 0.1, 1.0
@time tumor = birth_death_pushing(5000; b=b, d=d, μ=μ, ρ=Inf, dim=2)[:tumor] |> DataFrame

x = 0.3
tumor_reduced = reduced_μ(tumor, x)

using Plots, StatsBase

fig = tumor |> mutation_freqs |> df -> M(df.frequency; res=0.001, lab="")
Plots.plot!(collect(0:1:1/0.001), finv -> μ * b / (b - d) * (finv - 1), lab="")
tumor_reduced |> mutation_freqs |> df -> M!(fig, df.frequency; res=0.001, lab="")
Plots.plot!(collect(0:1:1/0.001), finv -> μ * x * b / (b - d) * (finv - 1), lab="")
