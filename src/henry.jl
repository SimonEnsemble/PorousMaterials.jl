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

    # get xtal density for conversion to per mass units (up here in case fails due to missing atoms in atomicmasses.csv)
    ρ = crystal_density(framework) # kg/m³

    # Bool's of whether to compute guest-host and/or guest-guest electrostatic energies
    #   there is no point in going through the computations if all charges are zero!
    const charged_system = charged(framework, verbose=verbose) & charged(molecule, verbose=verbose)

    # define Ewald summation params
    # pre-compute weights on k-vector contributions to long-rage interactions in
    #   Ewald summation for electrostatics
    #   allocate memory for exp^{i * n * k ⋅ r}
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sqrt(ljforcefield.cutoffradius_squared),
                                                          framework.box,
                                                          verbose=(verbose & charged_system),
                                                          ϵ=ewald_precision)

    # Σᵢ Eᵢe^(-βEᵢ)
    wtd_energy_sum = PotentialEnergy(0.0, 0.0)
    # Σᵢ e^(-βEᵢ)
    boltzmann_factor_sum = 0.0

    # conduct Monte Carlo insertions for less than 5 cores using Julia pmap function
    pmap(henryforloop, nb_insertions)

    # ⟨U⟩ = Σ Uᵢ e ^(βUᵢ) / [ ∑ e^(βUᵢ) ]
    average_energy = wtd_energy_sum / boltzmann_factor_sum

    result = Dict{String, Float64}()
    result["henry coefficient [mol/(m³-bar)]"] = boltzmann_factor_sum / (universal_gas_constant * temperature * nb_insertions)
    result["henry coefficient [mmol/(g-bar)]"] = result["henry coefficient [mol/(m³-bar)]"] / ρ
    result["henry coefficient [mol/(kg-Pa)]"] = result["henry coefficient [mmol/(g-bar)]"] / 100000.0
    result["⟨U⟩ (K)"] = average_energy.vdw + average_energy.coulomb
    result["⟨U, vdw⟩ (K)"] = average_energy.vdw
    result["⟨U, Coulomb⟩ (K)"] = average_energy.coulomb
    result["⟨U⟩ (kJ/mol)"] = result["⟨U⟩ (K)"] * 8.314 / 1000.0
    result["⟨U, vdw⟩ (kJ/mol)"] = result["⟨U, vdw⟩ (K)"] * 8.314 / 1000.0
    result["⟨U, Coulomb⟩ (kJ/mol)"] = result["⟨U, Coulomb⟩ (K)"] * 8.314 / 1000.0
    result["Qst (kJ/mol)"] = -result["⟨U⟩ (kJ/mol)"] + 8.314 * temperature / 1000.0

    if verbose
        print_with_color(:yellow, "\t----- results ----\n")
        for key in ["henry coefficient [mmol/(g-bar)]", "⟨U, vdw⟩ (kJ/mol)", "⟨U, Coulomb⟩ (kJ/mol)", "Qst (kJ/mol)"]
            @printf("\t %s = %f\n", key, result[key])
        end
    end
    return result
end

function henryforloop(framework::Framework, molecule::Molecule, temperature::Float64,
                           ljforcefield::LennardJonesForceField; nb_insertions::Int=1e6,
                           verbose::Bool=true, ewald_precision::Float64=1e-6)
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
end

function conduct_Widom_insertion(framework::Framework, w_molecule::Molecule, temperature::Float64,
                                ljforcefield::LennardJonesForceField,
                                verbose::Bool=true, ewald_precision::Float64=1e-6)

    repfactors = replication_factors(framework.box, ljforcefield)
    if verbose
        @printf("\tReplicating framework %d by %d by %d for short-range cutoff %.2f\n",
                repfactors..., sqrt(ljforcefield.cutoffradius_squared))
    end

    framework = replicate(framework, repfactors)
    ρ = crystal_density(framework) # kg/m³

    const charged_system = charged(framework, verbose=verbose) & charged(w_molecule, verbose=verbose)

    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sqrt(ljforcefield.cutoffradius_squared),
                                                          framework.box,
                                                          verbose=(verbose & charged_system),
                                                          ϵ=ewald_precision)

    x = framework.box.f_to_c * [rand(), rand(), rand()]
    translate_to!(w_molecule, x)
    if rotatable(w_molecule)
        rotate!(w_molecule)
    end

    energy = PotentialEnergy(0.0, 0.0)
    energy.vdw = vdw_energy(framework, w_molecule, ljforcefield)
    if charged_system
        energy.coulomb = electrostatic_potential_energy(framework, w_molecule, eparams,
                                                        kvecs, eikar, eikbr, eikcr)
    end

    boltzmann_factor = exp(-sum(energy) / temperature)
    return boltzmann_factor, energy
end
