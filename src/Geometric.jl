function unit_sphere_random_number()
    """
    Generate Cartesian coord on surface of unit sphere (faster)
    """
	ransq = 1.0
	ran1 = 0.0
	ran2 = 0.0
	while ransq >= 1.0
		ran1 = 2.0 * randn() - 1.0
		ran2 = 2.0 * randn() - 1.0
		ransq = ran1*ran1 + ran2*ran2
	end
	ranh = 2.0 * sqrt(1.0 - ransq)
	x = ran1 * ranh
	y = ran2 * ranh
	z = 1.0 - 2.0 * ransq
	return x, y, z
end

function unit_sphere_random_number2()
    """
    Generate Cartesian coord on surface of unit sphere (from cory's blog)
    """
	v = [0, 0, 0]
	while norm(v) < .0001
		x = randn()
		y = randn()
		z = randn()
		v = [x, y, z]
	end
	v = v / norm(v)
	return v[1], v[2], v[3]
end

function apply_boundary_conditions(distance, framework)
    """
    Apply periodic boundary conditions
    """
	x = zeros(Float64, 3)
	fractional = framework.cartesian_to_f_mtrx * [distance[1], distance[2], distance[3]]
	# apply boundary conditions
	x[1] = fractional[1] - round(fractional[1])
	x[2] = fractional[2] - round(fractional[2])
	x[3] = fractional[3] - round(fractional[3])
	cartesian = framework.f_to_cartesian_mtrx * [x[1], x[2], x[3]]
	return cartesian
end

function check_surface_area_overlap(probe, probe_size, f_xf, f_yf, f_zf, sigmas, framework)
    """
    Loop over all framework atoms, check for overlap (except self)
    """
	fAtom = zeros(Float64, 3)
	distance = zeros(Float64, 3)
	fAtom = [f_xf, f_yf, f_zf]
	for i = 1:framework.natoms
		fAtom_tmp = zeros(Float64, 3)
		fAtom_tmp = framework.fractional_coords[:, i]
		if fAtom != fAtom_tmp
			vdW_fAtom = sigmas[i]
			equilibrium_distance = 0.5 * (vdW_fAtom + probe_size)
			fAtom_xyz = framework.f_to_cartesian_mtrx * fAtom_tmp
			distance[1] = probe[1] - fAtom_xyz[1]
			distance[2] = probe[2] - fAtom_xyz[2]
			distance[3] = probe[3] - fAtom_xyz[3]
			dr = apply_boundary_conditions(distance, framework)
			rr = (dr[1] * dr[1]) + (dr[2] * dr[2]) + (dr[3] * dr[3])
			if rr < (equilibrium_distance * equilibrium_distance)
				return true
			end
		end
	end
	return false
end
'''
function vdw_energy(framework::Framework, molecule::Molecule,
                    ljforcefield::LennardJonesForceField, repfactors::Tuple{Int, Int, Int})
	energy = 0.0
    for ljsphere in molecule.ljspheres
        energy += vdw_energy(framework, ljsphere, ljforcefield, repfactors)
    end
	return energy
end
'''
function surface_area(framework::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, num_sampling=1000, probe_size=3.31)
	"""
	Calculates the accessible geometric surface area of an input structure
	"""

    epsilons, sigmas = _generate_epsilons_sigmas(framework, [forcefield])
	array_size = framework.natoms
	@printf("Size of probe: %f \n", probe_size)
	@printf("Unit Cell Volume: %f \n", framework.v_unitcell)
	@printf("Looping through atoms in %s to calculate surface area of %i atoms ... \n", structurename, framework.natoms)
	fractional_data = open("freq.txt", "w")
	SurfaceAreaAverage = 0.0
	for i = 1:array_size
        # print progress
        if i % (int(array_size / 10.0)) == 0
            @printf("\tPercent finished: %.1f\n", 100.0 * i / array_size)
        end
		total = 0.0
		counted = 0
		equilibrium_distance = 0.5 * (probe_size + sigmas[i])
		for j in 1:num_sampling
			total += 1
			probe = zeros(Float64, 3)
			x, y, z = unit_sphere_random_number2()
			fAtom_xyz = framework.f_to_cartesian_mtrx * framework.fractional_coords[:, i]
			probe[1] = fAtom_xyz[1] + x * equilibrium_distance
			probe[2] = fAtom_xyz[2] + y * equilibrium_distance
			probe[3] = fAtom_xyz[3] + z * equilibrium_distance
			# check for overlap with framework atoms
			overlap = check_surface_area_overlap(probe, probe_size, framework.fractional_coords[1, i], framework.fractional_coords[2, i], framework.fractional_coords[3, i], sigmas, framework)
			if overlap == false
				counted += 1
			end
		end
		tmp = (counted / total) * 4 * pi * (equilibrium_distance * equilibrium_distance)
		#@printf("%s %f %f %f %f \n", framework.atoms[i], counted, total, tmp, 1.0e4*tmp/framework.v_unitcell)
		SurfaceAreaAverage += tmp
		@printf(fractional_data, "%s %f %f %f \n", framework.atoms[i], counted / total, tmp, 1.0e4 * tmp / framework.v_unitcell)
	end
	@printf("Volumetric Surface Area (m2/cm3): %f", 1.0e4 * SurfaceAreaAverage / framework.v_unitcell)
    return 1.0e4 * SurfaceAreaAverage / framework.v_unitcell
end
