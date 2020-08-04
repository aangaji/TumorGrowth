using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise

@time using TumorGrowth

####### SIMULATION ######
#########################

@time (id, mut, t), tumor, mut_events = birth_death_pushing(2000; b=0.3, d=0.0, mu=0.1, dim=2)
tumor = tumor |> DataFrame

plotting_2d(tumor; annotate = false)
plotting_colored_mutations(tumor)

@show tumor
#########################
##### TUMOUR SAVING #####
#########################

using DataFrames, CSV

(id, mut, t), tumor, mut_events = birth_death_pushing(5000; b=0.3, d=0.0, mu=0.3, dim=2)
CSV.write("test_set_2d.csv", DataFrame(tumor), delim='\t')
b = data_import("test_set_2d.csv")
@show b

########################
####### PLOTTING #######
########################

b = data_import("test_set_3d.csv")
@time plotting(b)
@time plotting_colored_mutations(b)

b = data_import("test_set_2d.csv")
@time plotting_2d(b)
@time plotting_2d_colored_mutations(b)

########################
####### Sampling #######
########################

b = data_import("test_set_2d.csv")

fig = b |> plotting_colored_mutations
cross_section(b; x=10., width=3.) |> data -> plotting!(fig, data; color=:black)
cross_section(b; y=25., width=3.) |> data -> plotting!(fig, data; color=:black, path="test_set_2d")

fig = b |> plotting_colored_mutations
radial_sample(b; r=60., width=3.) |> data -> plotting!(fig, data; color=:black)

b = data_import("test_set_3d.csv")
scene = b |> plotting_colored_mutations
cross_section(b; x=10., width=3.) |> data -> plotting!(scene, data; color=:black, path="test_set_3d")
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
