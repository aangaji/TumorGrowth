using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise

@time using TumorGrowth

#########################
####### SIMULATION ######
#########################

@time (index, mut, t), tumor, mut_events = birth_death_pushing(10.; b=0.69, d=0.0, mu=0.3, dim=2, seed=1010)

birth_death_pushing!(tumor, mut_events, length(tumor)+1; b=0.69, d=0.0, mu=0.3, dim=2, t=t, cur_id=index, cur_mutation=mut, seed=1010)

b = tumor |> DataFrame

fig = plotting_colored_mutations(b, colormap=:tab20, color_set_size=20, shading=false, limits=limits)

# plotting_2d(tumor; annotate = false)
# plotting_2d_colored_mutations(tumor)
@show tumor

#########################
##### TUMOUR SAVING #####
#########################

using DataFrames, CSV

(index, mut, t), tumor, mut_events = birth_death_pushing(5000; b=0.3, d=0.0, mu=0.3, dim=2, seed=1234)
CSV.write("test_set_2d.csv", DataFrame(tumor), delim='\t')
b = data_import("test_set_2d.csv")
@show b

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
bulk(b; pos=(20.,20.), box=(20,30)) |> data -> plotting!(fig, data; color=:black)
punch(b; pos=(-20,-20), r=10) |> data -> plotting!(fig, data; color=:black)


b = data_import("test_5000_3d.csv")
scene = b |> plotting_colored_mutations
cross_section(b; x=10., width=3.) |> data -> plotting!(scene, data; color=:black, path="test_5000_3d.png")
cross_section(b; y=5., width=3.) |> data -> plotting!(scene, data; color=:black)

fig = cross_section(b; x=10., width=6.) |> plotting_colored_mutations
radial_sample(b; r=20., width=3.) |> data->cross_section(data; x=10., width=6.) |> data -> plotting!(fig, data; color=:black)


########################
####### Analysis #######
########################

b = data_import("test_set_3d.csv")

b |> plotting_colored_mutations
b.mutations |> allele_population
b.mutations |> mutation_freqs

########################
######## clones ########
########################

fig = plotting_colored_mutations(tumor, colormap=:tab20, color_set_size=20, shading=false)
clone(tumor, nothing) |> t -> plotting!(fig, t, color=:black)
clone(tumor, 1) |> t -> plotting!(fig, t, color=:black)

fig = plotting_colored_mutations(tumor, colormap=:tab20, color_set_size=20, shading=false)
clones(tumor)[1] |> t -> plotting!(fig, t, color=:black)

fig = plotting_colored_mutations(tumor, colormap=:tab20, color_set_size=20, shading=false)
clones(tumor)[10] |> clones .|> t -> plotting!(fig, t, color=:black)


#########################
###### time series ######
#########################

time_series = tumor_stepper(0.0:0.1:10.; b=0.69, d=0.0, mu=0.6, dim=2, seed = 1000)

show_clone!(scene, snapshot; mut) = plotting!(scene, clone(snapshot,mut); color=:black)

record_growth(time_series; path="test.gif", manipulate! = (sc, sn) -> show_clone!(sc, sn; mut=2))
