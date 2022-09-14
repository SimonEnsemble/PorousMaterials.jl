# add atom types to rc[:atomic_masses] for contextual LJ potentials
function append_atomic_masses()
    rc[:atomic_masses] = merge(rc[:atomic_masses], Dict(
        :N_in_N2 => 14.0067,
        :CH2 => 14.025,
        :CH3 => 15.035,
        :CH4 => 16.04,
        :C_b => 12.0107,
        :C_tol => 12.0107,
        :C_ac => 12.0107,
        :C_RCOO => 12.0107,
        :C_sp2 => 12.0107,
        :C_sp3 => 12.0107,
        :C_CO2 => 12.0107,
        :H_H2S => 1.00794,
        :H_b => 1.00794,
        :O_RCOO => 15.9994,
        :O_CO2 => 15.9994,
        :O_zeo => 15.9994,
        :S_H2S => 32.065,
        :Si_zeo => 28.0855,
        :ig => 1.0,
        :X => 1.0
    ))
end
