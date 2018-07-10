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
    # replicate the framework atoms so fractional coords are in [0, 1] spanning the simulation box
    framework = replicate(framework, repfactors)

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

    # conduct Monte Carlo insertions
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

        wtd_energy_sum += boltzmann_factor * energy
    end
    henry_coeff = boltzmann_factor_sum / (universal_gas_constant * temperature * nb_insertions) # mol/(m³-bar)
    average_energy = wtd_energy_sum / boltzmann_factor_sum

    if verbose
        println("henry coefficient (mol/(m³-bar)) = ", henry_coeff)
        println("⟨E, vdw⟩ (K) = ", average_energy.vdw)
        println("    (kJ/mol) = ", average_energy.vdw * 8.314 / 1000.0)
        println("⟨E, coulomb⟩ (K) = ", average_energy.coulomb)
        println("        (kJ/mol) = ", average_energy.coulomb * 8.314 / 1000.0)
    end
    return henry_coeff, average_energy
end
