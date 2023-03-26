using Plots
using Loess

include("extfns.jl")
includet("powerspectrum.jl")

function na_rm(x, y)
	newx = Float64[]
	newy = Float64[]

	for i in 1:length(y)
		if isnan(y[i]) != true
			push!(newx, x[i])
			push!(newy, y[i])
		end
	end
	return newx, newy
end

function graph_power_spectrum()
	k, P, Δ2 = power_spectrum()
	k, P = na_rm(k, P)
	
	model_kP = loess(k, P, span=0.1)
	us = range(k[10], k[length(k)], step=.1)
	vs = predict(model_kP, us)

	s = scatter(
		k, P, 
		xscale = :log10, 
		yscale = :log10,
		xlabel = "log₁₀k (h Mpc⁻¹)", 
		ylabel = "P(k) (h⁻³ Mpc³)",
		xlims = (10^-2, 10^1),
		ylims = (10^-2, 10^3),
		color = :grey,
		minorgrid = true,
		minorgridalpha = .3,
		legend = :false,
		size = (1000,1000),
		titlefontsize=16, 
		guidefontsize=16,
		tickfontsize=12,
		title = "1d spherically averaged power spectrum",
		markersize=4,
		markerstrokewidth=.1,
		markercolor = :grey,
		framestyle = :box,
	)
	l = plot!(
		us, 
		vs, 
		legend=false,
		linewidth=4,
		color = :black,
	)

	savefig("./media/spherical.png")
	plot!()
end


function graph_dimless_power_spectrum()
	k, P, Δ2 = power_spectrum()
	k, Δ2 = na_rm(k, Δ2)
	
	model = loess(k, Δ2, span=0.1)
	us = range(k[10], k[length(k)], step=.1)
	vs = predict(model, us)

	s = scatter(
		k, Δ2, 
		xscale = :log10, 
		yscale = :log10,
		xlabel = "log₁₀k (h Mpc⁻¹)", 
		ylabel = "Δ²(k)",
		xlims = (10^-2, 10^1),
		ylims = (10^-3, 10^-1),
		color = :grey,
		minorgrid = true,
		minorgridalpha = .3,
		legend = :false,
		size = (1000,1000),
		titlefontsize=16, 
		guidefontsize=16,
		tickfontsize=12,
		title = "1d dimensionless power spectrum",
		markersize=4,
		markerstrokewidth=.1,
		markercolor = :grey,
		framestyle = :box,
	)
	l = plot!(
		us, 
		vs, 
		legend=false,
		linewidth=4,
		color = :black,
	)

	savefig("./media/dimless.png")
	plot!()
end
