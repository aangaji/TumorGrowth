using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise

@time using TumorGrowth

#########################
####### SIMULATION ######
#########################


@time (index, mut, t), tumor = birth_death_pushing(3000; b=0.69, d=0.0, mu=0.3, ρc=2., dim=2, seed=1010)

bumor = deepcopy(tumor)

birth_death_pushing!(bumor, length(tumor)+1; b=0.69, d=0.0, mu=0.3, dim=2, t=t, cur_id=index, cur_mutation=mut, seed=1002)

birth_death_pushing(15.; b=0.69, d=0.0, mu=0.3, dim=3, seed=1010)

using Plots: palette

tumordf = tumor |> DataFrame
fig = plotting_colored_mutations(tumordf, colorpalette = palette(:tab20), shading=false, inline=false)


#########################
##### TUMOUR SAVING #####
#########################

using DataFrames, CSV

(index, mut, t), tumor = birth_death_pushing(2000; b=0.3, d=0.0, mu=0.3, ρc=1.8, dim=2, seed=1234)
CSV.write("test_2000_2d.csv", DataFrame(tumor), delim='\t')
b = data_import("test_2000_2d.csv")
@show b

b = data_import("test_2000_2d.csv")
@btime data_import("test_2000_2d.csv")

########################
####### PLOTTING #######
########################

b = data_import("test_2000_2d.csv")
@time plotting(b)
@time plotting_colored_mutations(b, shading = false)

b = data_import("test_5000_3d.csv")
@time plotting(b)
@time plotting_colored_mutations(b, shading=true)

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

b = data_import("test_5000_3d.csv")

b |> mutation_freqs

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

time_series = tumor_stepper(0.0:0.1:10.; b=0.69, d=0.0, mu=0.6, dim=2, seed = 1000)

function show_clone!(scene, snapshot; mut)
    plotting_colored_mutations!(scene, snapshot)
    plotting!(scene, clone(snapshot,mut); color=:black)
end

record_growth(time_series; path="test.gif", plot_func! = (sc, sn) -> show_clone!(sc, sn; mut=2))

#########################
###### density ~ b ######
#########################

d, ρc = 0.2, 2.
@time (index, mut, t), tumor = birth_death_pushing(5000; b=0.69, d=d, mu=0.3, ρc=ρc, dim=2, seed=1010)

using Plots
plot(-7.:0.01:7.,r-> TumorGrowth.w(r; σ=3.), fill=true)
plot!(-9.:0.01:9.,r-> TumorGrowth.w(r; σ=3.), xticks=-7:2:7, legend=:none)

tumordf = tumor |> DataFrame

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

@time tumor = birth_death_pushing(5000; b=0.69, d=0.1, mu=0.3, ρc=1.8, dim=2, seed=1010)[2] |> DataFrame

lattice, samples = multi_region_sampling(tumor; n = 100, cells_per_sample = 20)

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

lattice, samples = multi_region_sampling(tumor; n = 100, cells_per_sample = 20)
seq_results = filter!.(c->c.frequency > 0.0, mutation_freqs.(samples))

sampletumor = DataFrame(index = 1:length(samples), position = lattice, mutations = getproperty.(seq_results, :mutation), frequencies = getproperty.(seq_results, :frequency))

density = size(punch(tumor; pos = mean(tumor.position), r = 20.), 1)/(π*20^2)
scene = plotting_colored_mutations(sampletumor; markersize = sqrt(20/density/π), shading = true)
plotting_colored_mutations!(scene, tumor)

# save("multi_region_sequencing.png", scene)


##################################
## Reduce mutation artificially ##
##################################

b, d, μ = 0.69, 0.1, 1.0
@time tumor = birth_death_pushing(5000; b=b, d=d, mu=μ, ρc=Inf, dim=2)[2] |> DataFrame

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
