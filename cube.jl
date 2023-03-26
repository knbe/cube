using DelimitedFiles

const global data::Matrix{Float64} = readdlm("data/deltax_z8_50_cube_fix.dat", ',')

struct Cube
	N::Int64
	D::Int64
	Nd::Int64
	L::Float64
	V::Float64
	a::Float64
	δ::Array{Float64, 3}
end

function elle_arrive_de_la_mer()
	N = size(data)[1]
	D = 3
	Nd = round(Int, (N^(1/D)))
	L = 400.0
	V = L^D
	a = L/Nd	
	δ = zeros(Nd, Nd, Nd)
	
	for row in 1:N
		i = Int64(data[row,:][1] + 1)
		j = Int64(data[row,:][2] + 1)
		k = Int64(data[row,:][3] + 1)
		x = Float64(data[row,:][4])
		δ[i,j,k] = x
	end

	return Cube(N, D, Nd, L, V, a, δ)
end

const global cube = elle_arrive_de_la_mer()
