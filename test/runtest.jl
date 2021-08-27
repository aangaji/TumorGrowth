using Test

@testset "Simulation" begin

    simoutput = birth_death_pushing(100; b=1., d=0.0, μ=0.1, ρ=1., dim=2, seed=1010, showprogress=false)
    @test simoutput isa Dict
    @test map((:index,:mutation,:tumor,:time)) do key
        haskey(simoutput,key)
    end |> all
    @test try
        birth_death_pushing(100; b=1., d=1., μ=0.1, ρ=1., dim=2, seed=1010, showprogress=false)
    catch e
        e == ErrorException("Tumor died")
    end

    birth_death_pushing!(simoutput[:tumor], length(simoutput[:tumor])+1; b=0.69, d=0.0, μ=0.3, dim=2,
     t=simoutput[:time], cur_id=simoutput[:index], cur_mutation=simoutput[:mutation])

     tumor = simoutput[:tumor] |> DataFrame
     @test map(eachcol(tumor)) do col
         eltype(eltype(col)) <: Number
     end |> all

end

@testset "Plotting" begin
    tumor = birth_death_pushing(100; b=1., d=0.0, μ=0.1, ρ=1., dim=2, seed=1010, showprogress=false)[:tumor] |> DataFrame

    sc = plotting(tumor, shading=false, inline=false, size=(500,500))
    plotting_colored_mutations!(sc, tumor, colorpalette = TumorGrowth.palette(:tab20), show_warning = true)
    @test sc isa TumorGrowth.Scene
end
