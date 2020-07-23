using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise

@time using TumorGrowth

#########################
####### SIMULATION ######
#########################

@time a = birth_death_space_restricted(tumor_size=1000, b=0.3, d=0.1, mu=0.1)
@time a = birth_death_pushing(tumor_size=3000, b=0.3, d=0.0, mu=0.1)
plotting_2d(a[1]; annotate = false)
plotting(a[1])

@time a = birth_death_space_restricted_3D(tumor_size=500, b=0.3, d=0.0, mu=0.3)
@time a = birth_death_pushing_3D(tumor_size=3000, b=0.3, d=0.0, mu=0.3)
plotting(a[1])

@show a[1]

#########################
##### TUMOUR SAVING #####
#########################

using DataFrames, CSV

a = birth_death_pushing(tumor_size=5000, b=0.3, d=0.0, mu=0.3)
CSV.write("test_set_2d.csv", a[1], delim='\t')
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

#########################
##### Crosssections #####
#########################


b = data_import("test_set_2d.csv")

fig = b |> plotting_colored_mutations
cross_section(b; x=10., width=3.) |> data -> plotting!(fig, data; color=:black)
cross_section(b; y=25., width=3.) |> data -> plotting!(fig, data; color=:black, path="test_set_2d")

b = data_import("test_set_3d.csv")
scene = b |> plotting_colored_mutations
cross_section(b; x=10., width=3.) |> data -> plotting!(scene, data; color=:black, path="test_set_3d")
cross_section(b; y=5., width=3.) |> data -> plotting!(scene, data; color=:black)


########################
####### Analysis #######
########################
