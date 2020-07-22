using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise

@time using TumorGrowth

#########################
####### SIMULATION ######
#########################

@time A = birth_death_space_restricted(tumor_size=1000, b=0.3, d=0.1, mu=0.1)
@time A = birth_death_pushing(tumor_size=3000, b=0.3, d=0.0, mu=0.1)
@show A[1]
plotting_normal(A[1]; annotate = false)

@time a = birth_death_space_restricted_3D(tumor_size=500, b=0.3, d=0.0, mu=0.3)
@time a = birth_death_pushing_3D(tumor_size=3000, b=0.3, d=0.0, mu=0.3)
@show a
plotting_normal_3D(a[1])

plane_x = cross_section_3D_simulation(a[1])[1]
plotting_normal(plane_x)
plotting_colored_mutations(plane_x)



#########################
##### TUMOUR SAVING #####
#########################

using DataFrames, CSV

a = birth_death_pushing_3D(tumor_size=5000, b=0.3, d=0.0, mu=0.3)
CSV.write("pushing_test_set.csv", a[1], delim='\t')
b = data_import("pushing_test_set.csv")
@show b

########################
####### PLOTTING #######
########################
using Plots
@time plotting_normal_3D(b)
@time plotting_pushing_colored_mutations_3D(b)
plane_x = cross_section_3D_simulation(b)[1]
@time plotting_normal(plane_x)
@time plotting_colored_mutations(plane_x)


function circle(pos_x, pos_y, r)
    theta = LinRange(0, 2*pi, 100)
    return pos_x .+ r*sin.(theta), pos_y .+ r*cos.(theta)
end

plotly()
@time begin
    fig = Plots.plot()
    for p in [rand(2)*10 for _=1:1000]
        Plots.plot!(fig, circle(p..., 1), seriestype = [:shape,], lw = 0.5, c = palette(:tab10)[rand(1:10)], linecolor = :black, legend = false, fillalpha = 0.5, aspect_ratio = 1)
    end
    display(fig)
end

@time Plots.scatter(rand(1000)*10, rand(1000)*10, ms=10, alpha = 0.8, c = [palette(:tab10)[rand(1:10)] for _=1:100], legend = false, aspect_ratio = 1)
