using Pkg; Pkg.activate(pwd()); Pkg.instantiate()

using Plots, LsqFit, Statistics

using Revise

@time using TumorGrowth

b, d, ρ = log(2), 0.2, Inf
simoutput = birth_death_pushing(2000; b=b, d=d, μ=0.3, ρc=ρ, dim=2, seed=1010)

tumor = simoutput[:tumor]
times = simoutput[:times]
Ns = simoutput[:sizes]
t, N = last(times), last(Ns)

plot(times, log.(Ns))
plot!(times[end÷2:end], log.(Ns[end÷2:end]), legend=:none)


# N(t) = a*(t-t0)^n
# logN = n*log(t-t0) + log(a)
growth_model(logt; n, loga=0.0) = n.*logt .+ loga
exp_model(t; λ, n0=1, t0=0.0) = λ * (t.-t0)

fit_res = curve_fit( (t, p) -> growth_model(t; n = p[1], loga=p[2]),
          log.(times[end÷2:end]),
          log.(Ns[end÷2:end]),
          [1., 0.]
          )

fit_exp_res = curve_fit( (t, p) -> exp_model(t; λ = p[1], t0=p[2]),
          times[end÷2:end],
          log.(Ns[end÷2:end]),
          [1., 0.]
          )

plot()

confidence = (mean(abs.(fit_res.resid)) - mean(abs.(fit_exp_res.resid)))/(mean(abs.(fit_res.resid)) + mean(abs.(fit_exp_res.resid)))
if confidence < 0.
    plot!(log.(times[end÷2:end]), log.(Ns[end÷2:end]), lab="d=$d ρ=$ρ", legend=:topleft, ylabel=:size, xlabel=:time, alpha=0.3)
    plot!(log.(times), log.(Ns), lab=:none, alpha=0.3)
    plot!(log.(times[end÷2:end]), growth_model(log.(times[end÷2:end]); n=fit_res.param[1], loga=fit_res.param[2]), lab="fit", ylims=(0,log(N)), lw=2.)
    println("power law: ",fit_res.param, " confidence: ", -confidence)
else
    plot!(times[end÷2:end], log.(Ns[end÷2:end]), lab="d=$d ρ=$ρ", legend=:topleft, ylabel=:size, xlabel=:time, alpha=0.3)
    plot!(times, log.(Ns), lab=:none, alpha=0.3)
    plot!(times[end÷2:end], exp_model(times[end÷2:end]; λ = fit_exp_res.param[1], t0=fit_exp_res.param[2]), lab="fit", ylims=(0,log(N)), lw=2.)
    println("exponential: ",fit_exp_res.param, " confidence: ", confidence)
end
plot!()


conf2color(c::Float64) = c<0. ? cgrad(:blues)[-c] : cgrad(:reds)[c]

plot()

powers, powers_std, confidences = Float64[], Float64[], Float64[]
let d = 0.2
    ρs = [1.,1.25,1.5,1.75,2.]# ,2.25,2.5,3.,4.,5.,6.,7.,8.,9.,10.,12.]

    for ρ in ρs
        powersp = Float64[]
        confidencesp = Float64[]
        for _=1:10
            rerun = true
            simoutput = Dict()
            while rerun
                try simoutput = birth_death_pushing(2000; b=log(2), d=d, μ=0.3, ρc=ρ, dim=2)
                    rerun = false
                catch e
                end
            end
            times = simoutput[:times]
            Ns = simoutput[:sizes]
            fit_res = curve_fit( (t, p) -> growth_model(t; n = p[1], loga=p[2]),
                      log.(times[end÷2:end]),
                      log.(Ns[end÷2:end]),
                      [1., 0.]
                      )

            fit_exp_res = curve_fit( (t, p) -> exp_model(t; λ = p[1], t0=p[2]),
                      times[end÷2:end],
                      log.(Ns[end÷2:end]),
                      [1., 0.]
                      )
            confidence = (mean(abs.(fit_res.resid)) - mean(abs.(fit_exp_res.resid)))/(mean(abs.(fit_res.resid)) + mean(abs.(fit_exp_res.resid)))

            #power = c < 0. ? fit_res.param[1] : fit_exp_res.param[1]
            push!(powersp, fit_res.param[1])
            push!(confidencesp, confidence)
        end
        push!(powers, mean(powersp) )
        push!(powers_std, std(powersp))
        push!(confidences, mean(confidencesp) )
        println("ρ$ρ $(mean(powersp)) $(mean(confidencesp))")
    end
    scatter!(ρs, powers, yerror = powers_std,
        c= collect(confidences) .|> conf2color,
        marker=:x, lab="d=$d ρ=$ρ", legend=:bottomright)
end

# savefig("powerlaw_fit_power_truegrowthcurve.png")
