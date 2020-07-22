using Pkg; Pkg.activate(pwd())
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
@time a = birth_death_pushing_3D(tumor_size=300, b=0.3, d=0.0, mu=0.3)
@show a
plotting_3D_normal(a[1])

plane_x = cross_section_3D_simulation(a[1])[1]
plotting_normal(plane_x)



#########################
##### TUMOUR SAVING #####
#########################

using BSON, DataFrames, CSV

systime = time_ns()

a = birth_death_3D_pushing_fixed_cell_number(tumor_size=300, b=0.3, d=0.0, mu=0.3)
s = a[1]
BSON.@save "pushing_test_set.bson" s
CSV.write("pushing_test_set.csv", s, delim='\t')
b = BSON.load("pushing_test_set.bson")[:s]


########################
####### PLOTTING #######
########################
using BSON, Makie, Plots
b = BSON.load("pushing_test_set.bson")[:s]
@time plotting_3D_normal(b)
@time plotting_3D_pushing_colored_mutations(b)
plane_x = cross_section_3D_simulation(a[1])[1]
@time plotting_normal(plane_x)
@time plotting_colored_mutations(plane_x)

#### ONE PLOTTING FUNCTION SHOULD HANDLE DATA OF ANY DIM
Makie.meshscatter((randn(100), randn(100))..., ms=0.5, color=rand(1:30, 100), colormap=distinguishable_colors(68))


function circle(pos_x, pos_y, r)
    theta = LinRange(0, 2*pi, 100)
    return pos_x .+ r*sin.(theta), pos_y .+ r*cos.(theta)
end

plotly()
@time begin
    fig = plot()
    for p in [rand(2)*10 for _=1:1000]
        Plots.plot!(fig, circle(p..., 1), seriestype = [:shape,], lw = 0.5, c = palette(:tab10)[rand(1:10)], linecolor = :black, legend = false, fillalpha = 0.5, aspect_ratio = 1)
    end
    display(fig)
end

@time scatter(rand(1000)*10, rand(1000)*10, ms=10, alpha = 0.8, c = [palette(:tab10)[rand(1:10)] for _=1:100], legend = false, aspect_ratio = 1)
