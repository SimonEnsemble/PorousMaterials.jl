using Plots
using PorousMaterials

files = readdir("data/crystals")

cubelengths = []
i = 0

for file in files
	try
		frame = read_crystal_structure_file(file)
		strip_numbers_from_atom_labels!(frame)
		mol = Molecule(1, ["He"], zeros(3,1), [0.0])
		append!(cubelengths, maximum(frame.box.f_to_c * [1,1,1]))	
	catch e
		@printf("%s caused the following error:\n%s\n",file,e)
		i += 1
	end
end
@printf("Number of messed up files: %d\n",i)
histogram(cubelengths, nbins=100)
savefig("figure1.png")
