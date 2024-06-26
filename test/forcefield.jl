module Forcefield_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Forcefield Tests" begin
    ljforcefield = LJForceField("Dreiding"; r_cutoff=12.5, mixing_rules="Lorentz-Berthelot") # Dreiding
    # test reading of force field
    @test ljforcefield.pure_σ[:He] == 1.0
    @test ljforcefield.pure_ϵ[:Zn] == 12.0
    @test ljforcefield.σ²[:Zn][:He] == ((1.0 + 3.0) / 2)^2
    @test ljforcefield.ϵ[:He][:Zn] == sqrt(12.0 * 3.0)
    @test ljforcefield.ϵ[:He][:Zn] == ljforcefield.ϵ[:Zn][:He] # symmetry
    @test ljforcefield.σ²[:He][:Zn] == ljforcefield.σ²[:Zn][:He] # symmetry
    @test ljforcefield.r²_cutoff == 12.5^2

    # TODO test two other mixing rules
    # TODO test rep factors better

    # test calculation of replication factors required
    crystal = Crystal("test_structure.cif") # .cif
    strip_numbers_from_atom_labels!(crystal)
    rep_factors = replication_factors(crystal.box, ljforcefield)
    @test rep_factors == (25, 25, 25)

    # test check for force field coverage
    SBMOF1 = Crystal("SBMOF-1.cif")
    @test forcefield_coverage(SBMOF1.atoms, ljforcefield)
    @test forcefield_coverage(SBMOF1, ljforcefield)
    @test !forcefield_coverage(SBMOF1.atoms, LJForceField("bogus"))
end
end
