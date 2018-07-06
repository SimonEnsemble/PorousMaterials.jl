const universal_gas_constant = 8.3144598e-5 # m³-bar/(K-mol)

"""
TODO doc string
"""
function henry_coefficient(framework::Framework, molecule::Molecule, temperature::Float64, ljforcefield::LennardJonesForceField;
                           nb_insertions::Int=1e6, verbose::Bool=true, ewald_precision::Float64=1e-6)
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
    eparams, kvectors, eikar, eikbr, eikcr = setup_Ewald_sum(sqrt(ljforcefield.cutoffradius_squared), framework.box,
                        verbose=(verbose & charged_system),
                        ϵ=ewald_precision)

    henry_coeff = 0.0

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
        U_vdw = vdw_energy(framework, molecule, ljforcefield)
        U_coulomb = 0.0
        if charged_system
            U_coulomb = electrostatic_potential_energy(framework, molecule, eparams, kvectors,
                                                       eikar, eikbr, eikcr)
        end

        # increment henry coefficient
        henry_coeff += exp(-(U_vdw + U_coulomb) / temperature)
    end
    return henry_coeff / (universal_gas_constant * temperature * nb_insertions) # mol/(m³-bar)
end
