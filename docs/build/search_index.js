var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "functions.html#PorousMaterials.gcmc_simulation",
    "page": "Functions",
    "title": "PorousMaterials.gcmc_simulation",
    "category": "function",
    "text": "results, molecules = gcmc_simulation(framework, molecule, temperature, pressure,\n                                     ljforcefield; n_sample_cycles=5000,\n                                     n_burn_cycles=5000, sample_frequency=5,\n                                     verbose=false, molecules=Molecule[],\n                                     eos=:ideal)\n\nRuns a grand-canonical (μVT) Monte Carlo simulation of the adsorption of a molecule in a framework at a particular temperature and pressure using a Lennard Jones force field.\n\nA cycle is defined as max(20, number of adsorbates currently in the system) Markov chain proposals. Current Markov chain moves implemented are particle insertion/deletion and translation.\n\nArguments\n\nframework::Framework: the porous crystal in which we seek to simulate adsorption\ntemperature::Float64: temperature of bulk gas phase in equilibrium with adsorbed phase   in the porous material. units: Kelvin (K)\npressure::Float64: pressure of bulk gas phase in equilibrium with adsorbed phase in the   porous material. units: bar\nmolecule::Molecule: a template of the adsorbate molecule of which we seek to simulate   the adsorption\nljforcefield::LJForceField: the molecular model used to describe the   energetics of the adsorbate-adsorbate and adsorbate-host van der Waals interactions.\nn_burn_cycles::Int: number of cycles to allow the system to reach equilibrium before   sampling.\nn_sample_cycles::Int: number of cycles used for sampling\nsample_frequency::Int: during the sampling cycles, sample e.g. the number of adsorbed   gas molecules every this number of Markov proposals.\nverbose::Bool: whether or not to print off information during the simulation.\nmolecules::Array{Molecule, 1}: a starting configuration of molecules in the framework.\n\nNote that we assume these coordinates are Cartesian, i.e. corresponding to a unit box.\n\neos::Symbol: equation of state to use for calculation of fugacity from pressure. Default\n\nis ideal gas, where fugacity = pressure.\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.adsorption_isotherm",
    "page": "Functions",
    "title": "PorousMaterials.adsorption_isotherm",
    "category": "function",
    "text": "results = adsorption_isotherm(framework, molecule, temperature, pressures,\n                              ljforcefield; n_sample_cycles=100000,\n                              n_burn_cycles=10000, sample_frequency=25,\n                              verbose=false, molecules=Molecule[],\n                              ewald_precision=1e-6, eos=:ideal)\n\nRun a set of grand-canonical (μVT) Monte Carlo simulations in parallel. Arguments are the same as gcmc_simulation, as this is the function run in parallel behind the scenes. The only exception is that we pass an array of pressures. To give Julia access to multiple cores, run your script as julia -p 4 mysim.jl to allocate e.g. four cores. See Parallel Computing.\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.stepwise_adsorption_isotherm",
    "page": "Functions",
    "title": "PorousMaterials.stepwise_adsorption_isotherm",
    "category": "function",
    "text": "results = stepwise_adsorption_isotherm(framework, molecule, temperature, pressures,\n                              ljforcefield; n_sample_cycles=100000,\n                              n_burn_cycles=10000, sample_frequency=10,\n                              verbose=true, molecules=Molecule[],\n                              ewald_precision=1e-6, eos=:ideal)\n\nRun a set of grand-canonical (μVT) Monte Carlo simulations in series. Arguments are the same as gcmc_simulation, as this is the function run behind the scenes. An exception is that we pass an array of pressures. The adsorption isotherm is computed step- wise, where the ending configuration from the previous simulation (array of molecules) is passed into the next simulation as a starting point. The ordering of pressures is honored. By giving each simulation a good starting point, (if the next pressure does not differ significantly from the previous pressure), we can reduce the number of burn cycles required to reach equilibrium in the Monte Carlo simulation. Also see adsorption_isotherm which runs the μVT simulation at each pressure in parallel.\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.gcmc_result_savename",
    "page": "Functions",
    "title": "PorousMaterials.gcmc_result_savename",
    "category": "function",
    "text": "file_save_name = gcmc_result_savename(framework_name, molecule_species\n                                    ljforcefield_name, temperature, pressure,\n                                    n_burn_cycles, n_sample_cycles)\n\nDetermine the name of files saved during the GCMC simulation, be molecule positions or results. It uses many pieces of information from the simulation to ensure the file name accurately describes what it holds.\n\nArguments\n\nframework_name::AbstractString: The porous crystal being tested\nmolecule_species::Symbol: The molecule being tested inside the porous crystal\nljforcefield_name::AbstractString: The molecular model being used in this   simulation to describe intermolecular Van der Waals interactions\ntemperature::Float64: The temperature used in the simulation units: Kelvin (K)\npressure::Float64: The pressure used in the simulation units: bar\nn_burn_cycles::Int: The number of burn cycles used in this simulation\nn_sample_cycles::Int: The number of sample cycles used in this simulation\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.LJSphere",
    "page": "Functions",
    "title": "PorousMaterials.LJSphere",
    "category": "type",
    "text": "Data structure for a Lennard-Jones sphere, containing its species and position in  fractional coordinates.\n\nExample use\n\nljs = LJSphere(:C, [0.0, 0.0, 0.0])\n\nAttributes\n\nspecies::Symbol: atom species name, e.g. :C\nxf::Array{Float64, 1}: fractional coordinates, e.g. [1.0, 0.0, 4.0].\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.PtCharge",
    "page": "Functions",
    "title": "PorousMaterials.PtCharge",
    "category": "type",
    "text": "Point charge data structure indicates its charge and position in fractional coordinates.\n\nExample use\n\nptc = PtCharge(-0.2, [0.0, 0.0, 0.0])\n\nAttributes\n\nq::Float64: signed magnitude of charge (units: electrons), e.g. 1.0\nxf::Array{Float64, 1}: fractional coordinates, e.g. [1.0, 0.0, 4.0].\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.read_xyz",
    "page": "Functions",
    "title": "PorousMaterials.read_xyz",
    "category": "function",
    "text": "atoms, x = read_xyz(filename)\n\nReturn the list of atoms (Array{Symbol, 1}) and their Cartesian coordinates x::Array{Float64, 2} as stored in the .xyz file. x[:, k] will return Cartesian coords of the kth atom.\n\nArguments\n\nfilename::AbstractString: The filename of the .xyz file\n\nReturns\n\natoms::Array{Symbol, 1}: An array of atoms stored as symbols e.g. [:H, :H, :O] read\n\nfrom the .xyz file.\n\nx::Array{Float64, 2}: The Cartesian coordinates of the atoms. x[:, k] will return cartesian coordinates of the k-th atom\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.read_cpk_colors",
    "page": "Functions",
    "title": "PorousMaterials.read_cpk_colors",
    "category": "function",
    "text": "atom_colors = read_cpk_colors()\n\nRead in CPK color scheme for atoms. Return atom_colors::Dict{Symbol, Tuple{Int, Int, Int}} such that atom_colors[\":C\"] gives RGB code for carbon as a tuple, (144, 144, 144). https://en.wikipedia.org/wiki/CPK_coloring\n\nReturns\n\natom_colors::Dict{Symbol, Tuple{Int, Int, Int}}: A dictionary linking an element symbol to its\' corresponding CPK color in RGB\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.read_atomic_radii",
    "page": "Functions",
    "title": "PorousMaterials.read_atomic_radii",
    "category": "function",
    "text": "atomic_radii = read_atomic_radii()\n\nReturn atomic_radii::Dict{Symbol, Float64}, where atom_masses[\":C\"] gives the atomic radii of carbon (10.87 Angstrom).\n\nReturns\n\natomic_radii::Dict{Symbol, Float64}: A dictionary linking an element symbol to its\' corresponding atomic radius\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.write_to_xyz",
    "page": "Functions",
    "title": "PorousMaterials.write_to_xyz",
    "category": "function",
    "text": "write_to_xyz(atoms, x, filename; comment=\"\")\nwrite_to_xyz(molecules, box, filename; comment=\"\")\nwrite_to_xyz(framework, filename; comment=\"\")\n\nWrite a molecule, framework, or array of atoms & positions to an .xyz file.\n\nArguments\n\natoms::Array{Symbol, 1}: An array of atoms stored as symbols e.g. [:H, :H, :O]\nx::Array{Float64, 2}: The Cartesian coordinates of the atoms.\n\nx[:, k] contains Cartesian coordinates of the k-th atom\n\nmolecules::Array{Molecule, 1}: an array of molecules whose atoms to write to .xyz\nframework::Framework: a crystal structure whose atoms to write to .xyz\nfilename::AbstractString: The filename of the .xyz file. (\".xyz\" appended automatically\n\nif the extension is not provided.) (absolute path)\n\ncomment::AbstractString: comment if you\'d like to write to the file.\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.lennard_jones",
    "page": "Functions",
    "title": "PorousMaterials.lennard_jones",
    "category": "function",
    "text": "energy = lennard_jones(r², σ², ϵ)  (units: Kelvin)\n\nCalculate the lennard jones potential energy given the square of the radius r between two lennard-jones spheres. σ and ϵ are specific to interaction between two elements. Return the potential energy in units Kelvin (well, whatever the units of ϵ are).\n\nArguments\n\nr²::Float64: distance between two (pseudo)atoms in question squared (Angstrom²)\nσ²::Float64: sigma parameter in Lennard Jones potential squared (units: Angstrom²)\nϵ::Float64: epsilon parameter in Lennard Jones potential (units: Kelvin)\n\nReturns\n\nenergy::Float64: Lennard Jones potential energy\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.vdw_energy",
    "page": "Functions",
    "title": "PorousMaterials.vdw_energy",
    "category": "function",
    "text": "energy = vdw_energy(framework, molecule, ljforcefield)\n\nCalculates the van der Waals interaction energy between a molecule and a framework. Applies the nearest image convention to find the closest replicate of a specific atom.\n\nWARNING: it is assumed that the framework is replicated sufficiently such that the nearest image convention can be applied. See replicate.\n\nArguments\n\nframework::Framework: Crystal structure\nmolecule::Molecule: adsorbate (includes position/orientation/atoms)\nljforcefield::LJForceField: Lennard Jones force field\n\nReturns\n\nenergy::Float64: Van der Waals interaction energy\n\n\n\ngg_energy = vdw_energy(molecule_id, molecules, ljforcefield, simulation_box)\n\nCalculates van der Waals interaction energy of a single adsorbate molecules[molecule_id] with all of the other molecules in the system. Periodic boundary conditions are applied, using the nearest image convention.\n\nArguments\n\nmolecule_id::Int: Molecule ID used to determine which molecule in molecules we wish to calculate the guest-guest interactions\nmolecules::Array{Molecule, 1}: An array of Molecule data structures\nljforcefield::LJForceField: A Lennard Jones forcefield data structure describing the interactions between different atoms\nsimulation_box::Box: The simulation box for the computation.\n\nReturns\n\ngg_energy::Float64: The guest-guest interaction energy of molecules[molecule_id] with the other molecules in molecules\n\n\n\n"
},

{
    "location": "functions.html#PorousMaterials.vdw_energy_no_PBC",
    "page": "Functions",
    "title": "PorousMaterials.vdw_energy_no_PBC",
    "category": "function",
    "text": "Assumes unit cell box is a unit cube and no periodic boundary conditions are applied.\n\n\n\n"
},

{
    "location": "functions.html#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": "#Functions This page contains all of the functions exported by PorousMaterials. They are sorted by the .jl files they are found in.##Box.jl    Box\r\n    replicate\r\n    UnitCube##Crystal.jl    Framework\r\n    read_crystal_structure_file\r\n    remove_overlapping_atoms\r\n    strip_numbers_from_atom_labels!\r\n    write_unitcell_boundary_vtk\r\n    chemical_formula\r\n    molecular_weight\r\n    crystal_density\r\n    construct_box\r\n    replicate\r\n    read_atomic_masses\r\n    charged(::Framework, ::Bool)\r\n    write_cif\r\n    assign_charges##ElectrostaticsEnergetics.jl    electrostatic_potential\r\n    electrostatic_potential_energy\r\n    precompute_kvec_wts\r\n    setup_Ewald_sum\r\n    total##Energetics_Util.jl    PotentialEnergy\r\n    SystemPotentialEnergy##EOS.jl    PengRobinsonsGas\r\n    calculate_properties##Forcefield.jl    LJForcefield\r\n    replication_factors\r\n    check_forcefield_coverage##GCMC.jl    gcmc_simulation\r\n    adsorption_isotherm\r\n    stepwise_adsorption_isotherm\r\n    gcmc_result_savename##Grid.jl    Grid\r\n    apply_periodic_boundary_condition\r\n    write_cube\r\n    read_cube\r\n    energy_grid##Henry.jl    henry_coefficient\r\n    henry_result_savename##Matter.jl    LJSphere\r\n    PtCharge##MChelpers.jl    insert_molecule!\r\n    delete_molecule!\r\n    translate_molecule!\r\n    reinsert_molecule!\r\n    rotatable##Misc.jl    read_xyz\r\n    read_cpk_colors\r\n    read_atomic_radii\r\n    write_to_xyz##Molecules.jl    translate_to!\r\n    rotate!\r\n    rotation_matrix\r\n    rand_point_on_unit_sphere\r\n    charged(::Molecule, ::Bool)##NearestImage.jl    nearest_image!\r\n    nearest_r²\r\n    nearest_r##VdWEnergetics.jl    lennard_jones\r\n    vdw_energy\r\n    vdw_energy_no_PBC"
},

]}
