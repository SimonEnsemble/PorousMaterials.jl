using PorousMaterials

@printf("Starting script:\n")
files = readdir("data/crystals")

@printf("Iterating through %d crystal file structures\n",length(files))

for filename in files
	try
		frame = read_crystal_structure_file(filename)
	catch e
		@printf("%s caused the following error:\n%s\n",filename,e)
	end
end
