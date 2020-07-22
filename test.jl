using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise

@time using TumorGrowth

#########################
####### SIMULATION ######
#########################

@time a = birth_death_space_restricted(tumor_size=1000, b=0.3, d=0.1, mu=0.1)
@time a = birth_death_pushing(tumor_size=1000, b=0.3, d=0.0, mu=0.1)
plotting_normal(a[1]; annotate = false)

@time a = birth_death_space_restricted_3D(tumor_size=500, b=0.3, d=0.0, mu=0.3)
@time a = birth_death_pushing_3D(tumor_size=3000, b=0.3, d=0.0, mu=0.3)
plotting_normal_3D(a[1])

@show a[1]

#########################
##### TUMOUR SAVING #####
#########################

using DataFrames, CSV

a = birth_death_pushing_3D(tumor_size=5000, b=0.3, d=0.0, mu=0.3)
CSV.write("pushing_test_set_3d.csv", b, delim='\t')
b = data_import("pushing_test_set_3d.csv")
@show b

########################
####### PLOTTING #######
########################

b = data_import("pushing_test_set_3d.csv")
@time plotting_normal_3D(b)
@time plotting_pushing_colored_mutations_3D(b)

b = data_import("pushing_test_set_2d.csv")
@time plotting_normal(b)
@time plotting_colored_mutations(b)


#########################
##### Crosssections #####
#########################


b = data_import("pushing_test_set_2d.csv")

fig = b |> plotting_normal
cross_section(b; x=10., width=3.) |> data -> plotting_normal!(fig, data; color=:black)
cross_section(b; y=25., width=3.) |> data -> plotting_normal!(fig, data; color=:black)

b = data_import("pushing_test_set_3d.csv")

cross_section(b; x=10., width=3., reduce_dim=true)
