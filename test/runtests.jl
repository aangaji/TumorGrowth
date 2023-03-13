using Test, TumorGrowth

@testset "Simulation" begin

    simoutput = birth_death_pushing(100; b=1., d=0.0, mu=0.1, rho=1., dim=2, seed=1010, showprogress=false)
    @test simoutput isa Dict
    @test map((:index,:mutation,:tumor,:time)) do key
        haskey(simoutput,key)
    end |> all
    @test try
        birth_death_pushing(100; b=1., d=1., mu=0.5, rho=1., dim=2, seed=1010, showprogress=false)
    catch e
        e == ErrorException("Tumor died")
    end

    birth_death_pushing!(simoutput[:tumor], simoutput[:mutations], length(simoutput[:tumor])+1; dim=2,
     t=simoutput[:time], cur_id=simoutput[:index], cur_mutation=simoutput[:mutation])

     tumor = simoutput[:tumor] |> DataFrame
     @test map(eachcol(tumor)) do col
         eltype(eltype(col)) <: Number
     end |> all

     @test iszero(TumorGrowth.b_linear(1.1; b_max=1., rho_c=1.))
     reduced_mu!(tumor, 1/2)

     tumor = birth_death_pushing(2000; b=1., d=0.0, mu=0.3, rho=6., dim=3, seed=1234)[:tumor]
     TumorGrowth.CSV.write("temp.csv", DataFrame(tumor), delim=',')
     data_import("temp.csv", delim=",")
end

@testset "Sampling" begin
    tumor = data_import("temp.csv", delim=",")
    @test isa.( [
        cross_section(tumor; x=10., width=3.),
        radial_sample(tumor; r=30., width=3.),
        bulk(tumor; pos=(10.,0.,10.), box=(20,20,20)),
        punch(tumor; pos=(-10,-10,0), r=5)
        ], DataFrame) |> all
    @test isa.([
        mutation_freqs(tumor),
        stochastic_sequencing(tumor; readdepth=10)
        ], DataFrame) |> all

    plane = cross_section(tumor; x=0., width=3., reduce_dim=true)
    lattice, samples, sample_r = multi_region_sampling(plane; n = 20, cells_per_sample = 10)
    @test lattice isa Vector{Vector{Float64}} &&
        samples isa Vector{DataFrame}
    samples, sampletumor = multi_region_sequencing(plane;
         n = 20, cells_per_sample = 10, stochastic = true, readdepth=10)
    @test samples isa Vector{DataFrame} && sampletumor isa DataFrame
end

@testset "Plotting" begin
    tumor = data_import("temp.csv", delim=",")

    sc = plotting(tumor, shading=false, inline=false, size=(500,500))
    plotting_colored_mutations!(sc, tumor, colorpalette = TumorGrowth.palette(:tab20), show_warning = true, shading=false, path = "temp.png")

    plane = cross_section(tumor; x=0., width=3., reduce_dim=true)
    samples, sampletumor = multi_region_sequencing(plane; n = 20)
    plotting_colored_mutations(sampletumor; markersize=sampletumor.sample_r, show_axis=true)
    plotting_sampletumor_pies(sampletumor, first=10, colorpalette = TumorGrowth.palette(:tab10))
    plotting_sampletumor_pies(sampletumor, mutations = [1,10,55])

    @test sc isa TumorGrowth.Scene
    display(sc)
end

@testset "Timeseries" begin
    N, simparams... = (N=1000, b=1., d=0.0, mu=0.3, rho=6., dim=3, seed=1234)
    tumor = birth_death_pushing(N; simparams...)[:tumor] |> DataFrame

    time_series = tumor_stepper(range(0., last( tumor.t_birth), length=50); simparams...)

    @test last(time_series) == tumor

    record_growth(time_series; path="temp.gif",
            frames=1, shading=true,
            points_colors = t-> colors_by_mutations(t; colorpalette = TumorGrowth.palette(:tab20), show_warning=false)
            )
end
rm("temp.png")
rm("temp.csv")
rm("temp.gif")
