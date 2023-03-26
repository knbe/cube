using Statistics
using FFTW

include("cube.jl")
includet("old/coordinates.jl")

function power_spectrum()
	vraispacing = cube.L/cube.Nd
	samplingrate = (vraispacing/2π)^(-1)
	kf = fftfreq(cube.Nd, samplingrate)

	nksamples = 440
	ub = 10.0
	k = exp10.(range(log10(0.01), log10(ub), length=nksamples))

	# transform into fourier space
	lefou = (cube.a)^3 .* fft(cube.δ)
	pierrot = @. real(lefou * conj(lefou) / cube.V)

	#println(lefou)

	ccs = idkcubeset(nksamples)
	cube_kk_correlation!(cube, ccs, kf, k, ub)

#	println(typeof(kf))
#	println(typeof(ccs))
#	println(typeof(pierrot))
	
	# calculate the 2 pt correlation in fourier space
	P = two_pt_correlation(k, ccs, pierrot)
	Δ2 = @. P * (k^3 / 2π^2)

	return k, P, Δ2
end
