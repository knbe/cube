using Plots
using DelimitedFiles
#using AbstractFFTs
using FFTW

# 3d configuration space cube of data
const global data::Matrix{Float64} = readdlm("deltax.dat", ',')

nx::Int64 = 200
ny::Int64 = 200
nz::Int64 = 200

D::Float64 = 400.
N_ax::Int64 = 200

mutable struct Cube{T<:Int64}
	x::Vector{T}
	y::Vector{T}
	z::Vector{T}
end

function Cube()
	x = Int64[]
	y = Int64[]
	z = Int64[]
	return Cube(x,y,z)
end

function cubeset(num::Int64)
	cube = Vector{Cube}(undef, num)
	for i in 1:num
		cube[i] = Cube()
	end
	return cube
end

function get_nearest_arg(vec::Vector{Float64}, val::Float64)
	index = argmin(abs.(vec .- val))
end

function load()
	# volume of a voxel in Mpc
	dV = (D/N_ax)^3

	# DFT sample frequencies for each axis
	# i.e. k values
	samplefrequencies = 1 / (D / N_ax / 2π)
	kvals_x = fftfreq(N_ax, samplefrequencies)
	kvals_y = fftfreq(N_ax, samplefrequencies)
	kvals_z = fftfreq(N_ax, samplefrequencies)

	kvals_z = kvals_z[1:Int(N_ax/2)]

#	println(dV)
#	println(D/N_ax/2π)
#	println(samplefrequencies)
	#println(kvals_x)

	numpoints = 16
	Δk = 1.0e-2
	kmin = 3.5e-2
	kmax = 3.0
	kmag = exp10.(range(log10(kmin), log10(kmax), length=16))

	coordinates = cubeset(numpoints)

	for i in 1:N_ax
		for j in 1:N_ax
			for k in 1:length(kvals_z)
				kk = sqrt(
					kvals_x[i]^2 + 
					kvals_y[j]^2 + 
					kvals_z[k]^2
				)

				kk > kmax && continue

				index = get_nearest_arg(kmag, kk)

				push!(coordinates[index].x, i)
				push!(coordinates[index].y, j)
				push!(coordinates[index].z, k)
				
			end
		end
	end

	#println(data[1])
	cube = zeros(N_ax, N_ax, N_ax)

	#for row in 1:size(data)[1]
	for row in 1:size(data)[1]
		i = Int64(data[row,:][1] + 1)
		j = Int64(data[row,:][2] + 1)
		k = Int64(data[row,:][3] + 1)
		x = Float64(data[row,:][4])

		cube[i,j,k] = x
	end

	#println(typeof(cube))
	#cube::AbstractArray{Complex{Float64}} = cube
	#println(typeof(cube))
	cube_k = fft(cube) .* dV

	Pk = cube_k .* conj.(cube_k) ./ (D^3)
	Pk = real.(Pk)

	#println(coordinates[10].x)
	#println(size(Pk))
	P = zeros(length(kmag))
	for i in 1:length(kmag)
		P[i] = 0.0
		accumulator = Float64[]
		for j in 1:length(coordinates[i].x)
			push!(accumulator, 
				Pk[coordinates[i].x[j], 
				coordinates[i].y[j], 
				coordinates[i].z[j]]
			)
			#P[i] += mean(accumulator)
		end
		P[i] += mean(accumulator)
	end

	println(P)

	p1 = plot(
		kmag, P, 
		xscale = :log10, 
		yscale = :log10
	)
	p2 = plot(
		kmag, 
		P .*  (kmag.^3 ./ (2*π^2)), 
		xscale = :log10, 
		yscale = :log10
	)
	plot(
		p1, p2, 
		size=(1200,600)
	)

end
