module Forcefield_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Forcefield Tests" begin
    ljforcefield = LJForceField("Dreiding.csv", cutoffradius=12.5,
        mixing_rules="Lorentz-Berthelot") # Dreiding
    # test reading of force field
    @test ljforcefield.pure_σ[:He] == 1.0
    @test ljforcefield.pure_ϵ[:Zn] == 12.0
    @test ljforcefield.σ²[:Zn][:He] == ((1.0 + 3.0) / 2) ^ 2
    @test ljforcefield.ϵ[:He][:Zn] == sqrt(12.0 * 3.0)
    @test ljforcefield.ϵ[:He][:Zn] == ljforcefield.ϵ[:Zn][:He] # symmetry
    @test ljforcefield.σ²[:He][:Zn] == ljforcefield.σ²[:Zn][:He] # symmetry
    @test ljforcefield.cutoffradius_squared == 12.5 ^ 2

    # test calculation of replication factors required
    frame = Framework("test_structure.cif") # .cif
    strip_numbers_from_atom_labels!(frame)
    rep_factors = replication_factors(frame.box, ljforcefield)
    @test rep_factors == (25, 25, 25)

    # test check for force field coverage
    @test check_forcefield_coverage(Molecule("CO2"), ljforcefield)
    @test check_forcefield_coverage(Framework("SBMOF-1.cif"), ljforcefield)
    @test ! check_forcefield_coverage(Molecule("CO2"), LJForceField("bogus.csv"))
    @test ! check_forcefield_coverage(Framework("SBMOF-1.cif"), LJForceField("bogus.csv"))
end
end
