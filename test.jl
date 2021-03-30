#=
        This test file demonstrates most functionality of the package.
        Load Revise if you want to make changes to the source files.
=#

using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise

@time using TumorGrowth

#########################
######## SANDBOX ########
#########################



#########################
####### SIMULATION ######
#########################

@time simoutput = birth_death_pushing(5000; b=0.69, d=0.0, μ=0.1, ρc=Inf, dim=2)#, seed=1010)
tumor = simoutput[:tumor]

@time simoutput = nonspatial(1000000; b=0.69, d=0.0, μ=0.3, dim=2)

bumor = deepcopy(tumor)

birth_death_pushing!(bumor, length(tumor)+1; b=0.69, d=0.0, μ=0.3, dim=2,
 t=simoutput[:time], cur_id=simoutput[:index], cur_mutation=simoutput[:mutation])

birth_death_pushing(15.; b=0.69, d=0.0, μ=0.3, dim=3, seed=1010)

using Plots: palette

tumordf = tumor |> DataFrame
fig = plotting_colored_mutations(tumordf, colorpalette = palette(:tab20), shading=false, inline=false)


#########################
##### TUMOUR SAVING #####
#########################

using DataFrames, CSV

tumor = birth_death_pushing(5000; b=0.3, d=0.0, μ=0.3, ρc=6., dim=3, seed=1234)[:tumor]
CSV.write("test_5000_3d.csv", DataFrame(tumor), delim='\t')
b = data_import("test_2000_2d.csv")
@show b

b = data_import("test_2000_2d.csv")
@btime data_import("test_2000_2d.csv")

########################
####### PLOTTING #######
########################

b = data_import("test_2000_2d.csv")
@time plotting(b)
@time plotting_colored_mutations(b, shading = false, path="test_2000_2d.png")

b = data_import("test_5000_3d.csv")
@time plotting(b)
@time plotting_colored_mutations(b, shading=true, path="test_5000_3d.png")

########################
####### Sampling #######
########################

b = data_import("test_2000_2d.csv")

fig = b |> plotting_colored_mutations
cross_section(b; x=10., width=3.) |> data -> plotting!(fig, data; color=:black)
cross_section(b; y=25., width=3.) |> data -> plotting!(fig, data; color=:black)
radial_sample(b; r=30., width=3.) |> data -> plotting!(fig, data; color=:black, path="test_2000_2d.png")

fig = b |> plotting_colored_mutations
bulk(b; pos=(25.,0.), box=(20,30)) |> data -> plotting!(fig, data; color=:black)
punch(b; pos=(-20,-20), r=10) |> data -> plotting!(fig, data; color=:black)


b = data_import("test_5000_3d.csv")
scene = b |> plotting_colored_mutations
cross_section(b; x=10., width=3.) |> data -> plotting!(scene, data; color=:black, path="test_5000_3d.png")
cross_section(b; y=5., width=3.) |> data -> plotting!(scene, data; color=:black)

display(scene)

fig = cross_section(b; x=10., width=6.) |> plotting_colored_mutations
radial_sample(b; r=20., width=3.) |> data->cross_section(data; x=10., width=6.) |> data -> plotting!(fig, data; color=:black)

########################
####### Analysis #######
########################

tumor = data_import("test_5000_3d.csv")

tumor |> mutation_freqs
tumor |> stochastic_sequencing

function M!(fig, f; res=0.0, nBins = 100, lab=:none, color=:auto, normalize=false)
        hist = fit(Histogram, 1 ./ f[res.<f], nbins=nBins, closed=:left)
        M = [ sum(hist.weights[1:n]) for n=1:length(hist.weights)]
        M = normalize ? M ./ last(M) : M
        plot!(fig, midpoints(hist.edges[1]), M,
                marker=(:o), ms=1.5, lab=lab, c=color, xlabel="1/f", ylabel="M(f)", legend=:bottomright)
end
M(f; res=0.0, nBins = 100, lab=:none, color=:auto, normalize=false) = M!(plot(), f; res=res, nBins = nBins, lab=lab, color=color, normalize=normalize)

using Plots, StatsBase, Statistics
M(mutation_freqs(tumor).frequency; res=0.0005)
b, d, μ=0.69, 0.00, 00.3
plot!(1:1000, n->μ*b/(b-d)*n )

f = filter(m -> m.reads > 1, stochastic_sequencing(tumor, readdepth=1000)).frequency
M(stochastic_sequencing(tumor, readdepth=2000).frequency; res=0.001)

########################
######## clones ########
########################

fig = plotting_colored_mutations(b; colorpalette=palette(:tab20), shading=false)
clone(b, nothing) |> t -> plotting!(fig, t, color=:black)
clone(b, 1) |> t -> plotting!(fig, t, color=:black)

fig = plotting_colored_mutations(b; colorpalette=palette(:tab20), shading=false)
clones(b)[1] |> t -> plotting!(fig, t, color=:black)

fig = plotting_colored_mutations(b; colorpalette=palette(:tab20), shading=false)
clones(b)[10] |> clones .|> t -> plotting!(fig, t, color=:black)


#########################
###### time series ######
#########################

using DataFrames, Makie
using Plots: palette, distinguishable_colors

time_series = tumor_stepper(0.0:0.1:40.; b=0.69, d=0.0, μ=0.3, ρc=1.7, dim=2, seed = 1000)

record_growth(time_series; path="test.gif",
        frames=1, shading=false,
        points_colors = t-> colors_by_mutations(t; colorpalette = distinguishable_colors(64))
        )


#########################
###### density ~ b ######
#########################

using Plots

TumorGrowth.b_linear()
TumorGrowth.b_hill(6)
TumorGrowth.b_curve(1.; bup=1., ρc=1.1)

d, ρc = 0.0, 1.
@time tumor = birth_death_pushing(5000; b=0.69, d=d, μ=0.3, ρc=ρc, dim=2)[:tumor]
tumordf = tumor |> DataFrame

# power law / surface growth
plot(log.(tumordf.t_birth), log.(1:length(tumor)), legend=:none )
plot!(log.(tumordf.t_birth), t -> 2*t, xlims=(0., 4.))

# exponential / bulk
plot!(tumordf.t_birth, log.(1:5000), legend=:none )
plot!(tumordf.t_birth, t -> (log(2) - d)*t )

plot(-7.:0.01:7.,r-> TumorGrowth.w(r; σ=3.), fill=true)
plot!(-9.:0.01:9.,r-> TumorGrowth.w(r; σ=3.), xticks=-7:2:7, legend=:none)

TumorGrowth.b_linear()
TumorGrowth.b_hill(6)
plot(0.:0.01:9., ρ-> TumorGrowth.b_curve(ρ; bup=1., ρc=1.), legend=:none)
vline!([1.])

plotting_colored_mutations(tumordf, colorpalette = palette(:tab20), shading=false, inline=false)

histogram(tumordf.b, alpha=0.3, lw=0.2, lab="b")
vline!([sum(tumordf.b)/length(tumor)], lab="mean", c=:blue)
vline!([d], lab="d", c=:black)

plotting(tumordf; color = tumordf.b./maximum(tumordf.b), colormap = cgrad(:reds))


#############################
### mulit region sampling ###
#############################

using Makie: meshscatter, meshscatter!, save, Point
using Statistics

@time tumor = birth_death_pushing(5000; b=0.69, d=0.1, μ=0.3, ρc=1.8, dim=2, seed=1010)[:tumor] |> DataFrame

lattice, samples = multi_region_sampling(tumor; n = 100, cells_per_sample = 20)

samples, sampletumor = multi_region_sequencing(tumor; n = 100, cells_per_sample = 20, stochastic = true)

mean(size.(samples,1))

scene = plotting_colored_mutations(tumor)
for sample in samples
    plotting!(scene, sample; color=:black)
end
display(scene)

sample_r = 2*sqrt(20/π)
cm = mean(tumor.position)
density = size(punch(tumor; pos=cm, r=10.), 1)/(π*10. ^2)
R = sqrt(size(tumor,1)/(π*density))
scene = plotting_colored_mutations(tumor)
meshscatter!([Point{2}.(lattice)...], markersize = sample_r, scale_plot = false)
meshscatter!(0:0.2:360 .|> ϕ-> Point{2}(R .*(cosd(ϕ),sind(ϕ)) .+ cm), markersize = 1. )

# save("multi_region_sampling.png", scene)

density = size(punch(tumor; pos = mean(tumor.position), r = 20.), 1)/(π*20^2)
scene = plotting_colored_mutations(sampletumor; markersize = sqrt(20/density/π), shading = true)
plotting_colored_mutations!(scene, tumor)

# save("multi_region_sequencing.png", scene)

##################################
## Reduce mutation artificially ##
##################################

b, d, μ = 0.69, 0.1, 1.0
@time tumor = birth_death_pushing(5000; b=b, d=d, μ=μ, ρc=Inf, dim=2)[:tumor] |> DataFrame

x = 0.1
tumor_reduced = reduced_μ(tumor, x)

function M!(fig, f; res=0.0, nBins = 100, lab=:none, color=:auto)
        hist = fit(Histogram, 1 ./ f[res.<f], nbins=nBins, closed=:left)
        M = [ sum(hist.weights[1:n]) for n=1:length(hist.weights)]
        plot!(fig, midpoints(hist.edges[1]), M,
                marker=(:o), ms=1.5, lab=lab, c=color, xlabel="1/f", ylabel="M(f)", legend=:bottomright)
end
M(f; res=0.0, nBins = 100, lab=:none, color=:auto) = M!(plot(), f; res=res, nBins = nBins, lab=lab, color=color)

using Plots: plot, plot!
using StatsBase

fig = tumor |> mutation_freqs |> df -> M(df.frequency; res = 0.001)
plot!(collect(0:1:1/0.001), finv -> μ*b/(b-d) *(finv - 1), lab=:none)
tumor_reduced |> mutation_freqs |> df -> M!(fig, df.frequency; res = 0.001)
plot!(collect(0:1:1/0.001), finv -> μ*x*b/(b-d) *(finv - 1), lab=:none)

plotting_colored_mutations(tumor)

plotting_colored_mutations(tumor_reduced)
