using ProgressMeter

const universal_gas_constant = 8.3144598e-5 # m³-bar/(K-mol)

"""
TODO doc string

TODO break into blocks and report standard deviations
TODO parallelize
"""
function henry_coefficient(framework::Framework, molecule::Molecule, temperature::Float64,
                           ljforcefield::LennardJonesForceField; nb_insertions::Int=1e6,
                           verbose::Bool=true, ewald_precision::Float64=1e-6)
    if verbose
        print("Simulating Henry coefficient of ")
        print_with_color(:green, molecule.species)
        print(" in ")
        print_with_color(:green, framework.name)
        print(" at ")
        print_with_color(:green, temperature)
        print(" K, using ")
        print_with_color(:green, ljforcefield.name)
        print(" force field with ")
        print_with_color(:green, nb_insertions)
        println(" insertions.")
    end

    # replication factors for applying nearest image convention for short-range interactions
    repfactors = replication_factors(framework.box, ljforcefield)
    if verbose
        @printf("\tReplicating framework %d by %d by %d for short-range cutoff %.2f\n",
                repfactors..., sqrt(ljforcefield.cutoffradius_squared))
    end
    # replicate the framework atoms so fractional coords are in [0, 1] spanning the simulation box
    framework = replicate(framework, repfactors)
    
    # Bool's of whether to compute guest-host and/or guest-guest electrostatic energies
    #   there is no point in going through the computations if all charges are zero!
    const charged_system = charged(framework, verbose=verbose) & charged(molecule, verbose=verbose)

    # get xtal density for conversion to per mass units (up here in case fails due to missing atoms in atomicmasses.csv)
    ρ = crystal_density(framework) # kg/m³
    
    # partition total insertions among blocks.
    if nprocs() > N_BLOCKS
        error("Use $N_BLOCKS cores or less for Henry coefficient calculations to match the number of blocks")
    end
    nb_insertions_per_block = ceil(Int, nb_insertions / N_BLOCKS)

    # conduct Monte Carlo insertions for less than 5 cores using Julia pmap function
    # set up function to take a tuple of arguments, the number of insertions to 
    # perform and the molecule to move around/rotate. each core needs a different
    # molecule because it will change its attributes in the simulation
    # x = (nb_insertions, molecule) for that core
    henry_loop(x::Tuple{Int, Molecule}) = _conduct_Widom_insertions(framework, x[2], 
                                            temperature, ljforcefield, x[1], 
                                            charged_system, ewald_precision, verbose)
    
    # parallelize insertions across the cores
    res = pmap(henry_loop, [(nb_insertions_per_block, deepcopy(molecule)) for b = 1:N_BLOCKS])

    # unpack the boltzmann factor sum and weighted energy sum from each block
    boltzmann_factor_sums = [res[b][1] for b = 1:N_BLOCKS] # Σᵢ e^(-βEᵢ) for that core
    wtd_energy_sums = [res[b][2] for b = 1:N_BLOCKS] # Σᵢ Eᵢe^(-βEᵢ) for that core
    
    # compute block ⟨U⟩, Kₕ
    #   ⟨U⟩ = Σ Uᵢ e ^(βUᵢ) / [ ∑ e^(βUᵢ) ]
    #   Kₕ = β Σ e ^(βUᵢ) / nb_insertions_per_block
    #    (these N_BLOCKS-long arrays)
    average_energies = wtd_energy_sums ./ boltzmann_factor_sums # K
    henry_coefficients = boltzmann_factor_sums / (universal_gas_constant * temperature * nb_insertions_per_block) # mol/(m³-bar)
    if verbose
        for b = 1:N_BLOCKS
            print_with_color(:yellow, @sprintf("\tBlock  %d/%d statistics:\n", b, N_BLOCKS))
            println("\tHenry coeff. [mmol/(g-bar)]: ", henry_coefficients[b] / ρ)
            println("\t⟨U, vdw⟩ (K): ", average_energies[b].vdw)
            println("\t⟨U, Coulomb⟩ (K): ", average_energies[b].coulomb)
        end
    end

    # compute error estimates
    err_kh = 2.0 * std(henry_coefficients) / sqrt(N_BLOCKS)
    err_energy = PotentialEnergy()
    err_energy.vdw = 2.0 * std([average_energies[b].vdw for b = 1:N_BLOCKS]) / sqrt(N_BLOCKS)
    err_energy.coulomb = 2.0 * std([average_energies[b].coulomb for b = 1:N_BLOCKS]) / sqrt(N_BLOCKS)

    result = Dict{String, Float64}()
    result["henry coefficient [mol/(m³-bar)]"] = mean(henry_coefficients)
    result["henry coefficient [mmol/(g-bar)]"] = result["henry coefficient [mol/(m³-bar)]"] / ρ
    result["err henry coefficient [mmol/(g-bar)]"] = err_kh / ρ
    result["henry coefficient [mol/(kg-Pa)]"] = result["henry coefficient [mmol/(g-bar)]"] / 100000.0
    
    # note assumes same # insertions per core.
    result["⟨U, vdw⟩ (K)"] = mean([average_energies[b].vdw for b = 1:N_BLOCKS])
    result["⟨U, Coulomb⟩ (K)"] = mean([average_energies[b].coulomb for b = 1:N_BLOCKS])
    result["⟨U⟩ (K)"] = result["⟨U, vdw⟩ (K)"] + result["⟨U, Coulomb⟩ (K)"]

    result["⟨U⟩ (kJ/mol)"] = result["⟨U⟩ (K)"] * 8.314 / 1000.0
    result["⟨U, vdw⟩ (kJ/mol)"] = result["⟨U, vdw⟩ (K)"] * 8.314 / 1000.0
    result["err ⟨U, vdw⟩ (kJ/mol)"] = err_energy.vdw * 8.314 / 1000.0
    result["⟨U, Coulomb⟩ (kJ/mol)"] = result["⟨U, Coulomb⟩ (K)"] * 8.314 / 1000.0
    result["err ⟨U, Coulomb⟩ (kJ/mol)"] = err_energy.coulomb * 8.314 / 1000.0
    result["Qst (kJ/mol)"] = -result["⟨U⟩ (kJ/mol)"] + 8.314 * temperature / 1000.0
    result["err Qst (kJ/mol)"] = sum(err_energy) * 8.314 / 1000.0

    if verbose
        print_with_color(:green, "\t----- final results ----\n")
        for key in ["henry coefficient [mmol/(g-bar)]", "⟨U, vdw⟩ (kJ/mol)", "⟨U, Coulomb⟩ (kJ/mol)", "Qst (kJ/mol)"]
            @printf("\t%s = %f +/- %f\n", key, result[key], result["err " * key])
        end
    end
    return result
end

# assumed framework is already replicated sufficiently for short-range interactions
# to facilitate parallelization
function _conduct_Widom_insertions(framework::Framework, molecule::Molecule, 
                                   temperature::Float64, ljforcefield::LennardJonesForceField,
                                   nb_insertions::Int, charged_system::Bool,
                                   ewald_precision::Float64, verbose::Bool)
    # define Ewald summation params; this must be done on each core since they cannot share eikar for race condition possibility
    # pre-compute weights on k-vector contributions to long-rage interactions in
    #   Ewald summation for electrostatics
    #   allocate memory for exp^{i * n * k ⋅ r}
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sqrt(ljforcefield.cutoffradius_squared),
                                                          framework.box,
                                                          verbose=(verbose & (myid() == 1)  & charged_system),
                                                          ϵ=ewald_precision)
    # to be Σᵢ Eᵢe^(-βEᵢ)
    wtd_energy_sum = PotentialEnergy(0.0, 0.0)
    # to be Σᵢ e^(-βEᵢ)
    boltzmann_factor_sum = 0.0

    for i = 1:nb_insertions
        # determine uniform random center of mass
        x = framework.box.f_to_c * [rand(), rand(), rand()]
        # translate molecule to the new center of mass
        translate_to!(molecule, x)
        # rotate randomly
        if rotatable(molecule)
            rotate!(molecule)
        end

        # calculate potential energy of molecule at that position and orientation
        energy = PotentialEnergy(0.0, 0.0)
        energy.vdw = vdw_energy(framework, molecule, ljforcefield)
        if charged_system
            energy.coulomb = electrostatic_potential_energy(framework, molecule, eparams,
                                                            kvecs, eikar, eikbr, eikcr)
        end

        # calculate Boltzmann factor e^(-βE)
        boltzmann_factor = exp(-sum(energy) / temperature)

        boltzmann_factor_sum += boltzmann_factor

        # to avoid NaN; contribution to wtd energy sum when energy = Inf is zero.
        if (! isinf(energy.vdw)) && (! isinf(energy.coulomb))
            wtd_energy_sum += boltzmann_factor * energy
        end
        # else add 0.0 b/c lim E --> ∞ E exp(-E) is zero.
    end
    return boltzmann_factor_sum, wtd_energy_sum
end
