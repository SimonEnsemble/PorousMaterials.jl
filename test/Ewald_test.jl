using PorousMaterials
using Base.Test
using BenchmarkTools, Compat

@printf("------------------------------\nTesting Ewald.jl\n\n")
const framework = read_crystal_structure_file("NU-1000_Greg.cif")

k_rep_factors = (11, 11, 9)
α = 0.265058
sr_cutoff = 12.5
rep_factors = replication_factors(framework, sr_cutoff)
sim_box = replicate_box(framework.box, rep_factors)
const kvectors = compute_kvectors(sim_box, k_rep_factors, α)

q_test = 0.8096

@testset "Ewald summation Tests" begin
    x = [9.535619863743, 20.685576379935, 0.127344239990]
    ϕ = electrostatic_potential(framework, x, sim_box, rep_factors, sr_cutoff, kvectors, α)
    @test isapprox(ϕ * q_test, 111373.38, atol=2.5)

    x = [4.269654927228, 23.137319129548, 28.352847101096]
    ϕ = electrostatic_potential(framework, x, sim_box, rep_factors, sr_cutoff, kvectors, α)
    @test isapprox(ϕ * q_test, -531.0, atol=0.5)

    x = [-0.047382031804, 7.209555961450, 5.158180463556]
    ϕ = electrostatic_potential(framework, x, sim_box, rep_factors, sr_cutoff, kvectors, α)
    @test isapprox(ϕ * q_test, -2676.8230141, atol=0.5)
end

# timing
x = [-0.047382031804, 7.209555961450, 5.158180463556]
@btime ϕ = electrostatic_potential(framework, x, sim_box, rep_factors, sr_cutoff, kvectors, α)
 # ϕ = ϕ_sr(framework, x, rep_factors, sr_cutoff, α)
 # @btime ϕ_sr(framework, x, rep_factors, sr_cutoff, α)
 # 
 # ϕ = ϕ_lr(framework, x, sim_box, rep_factors, kvectors, α)
 # @btime ϕ_lr(framework, x, sim_box, rep_factors, kvectors, α)
