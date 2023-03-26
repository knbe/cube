using GLMakie
using Colors

include("data.jl")

function show_cube(ρ)
	grid = [Point3f(i,j,k) for i in 1:200 for j in 1:200 for k in 1:200]
	cgrid = [RGBA(4 * ρ[i,j,k], 3.5 * ρ[i,j,k], 3.5 * ρ[i,j,k], 0.02) 
			for i in 1:200 for j in 1:200 for k in 1:200]

	ρgrid = [ρ[i,j,k] for i in 1:200 for j in 1:200 for k in 1:200]

	perspectiveness=0.4
	aspect = :data

	fontsize_theme = Theme(fontsize=30)
	set_theme!(fontsize_theme)

	fig = Figure(resolution=(1200,1200))
	ax = Axis3(fig[1,1]; aspect, perspectiveness, azimuth=.72)
	meshscatter!(ax, grid; color=cgrid, markersize=.5, shading=false)
	#meshscatter!(ax, grid; color=ρgrid, colormap=:plasma, markersize=.5, shading=false, alpha = 0.02)
	#meshscatter!(ax, grid; color=cgrad(:plasma, ρgrid, alpha = 0.05), shading=false)
	#limits!(ax, 200, 200, 200, 200,200,200,200,200,200)
	display(fig)
end
