#using GLVisualize
#using GLWindow
@printf("This shouldn't run")
function visualize_landscape(energy::Array{Float64, 3}, isovalue::Float64)

	if !isdefined(:runtests)
		window = glscreen()
		timesignal = bounce(linspace(0f0, 1f0, 360))
	end

	max = maximum(abs.(energy))
	min = minimum(energy)
	energy2 = (energy) ./ (max .- min)

	description = """
	Some text
	"""

	vol = visualize(energy2, :iso, isovalue=isovalue)
	_view(vol, window)

	if !isdefined(:runtests)
		renderloop(window)
	end
end
