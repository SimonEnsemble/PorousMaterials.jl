using PorousMaterials

f = read_crystal_structure_file("NU-1000_Greg.cif")

k_rep_factors = (11, 11, 9)
α = 0.265058
sr_cutoff = 12.5
rep_factors = replication_factors(f, sr_cutoff)
sim_box = replicate_box(f.box, rep_factors)

q_test = 0.8096

x = [9.535619863743, 20.685576379935, 0.127344239990]
ϕ = electrostatic_potential(f, x, sim_box, rep_factors, sr_cutoff, k_rep_factors, α)
println("ϕ = ", ϕ  * q_test)
println("should be 111373.38769648 K")

x = [4.269654927228, 23.137319129548, 28.352847101096]
ϕ = electrostatic_potential(f, x, sim_box, rep_factors, sr_cutoff, k_rep_factors, α)
println("ϕ = ", ϕ  * q_test)
println("should be -531 K")

x = [-0.047382031804, 7.209555961450, 5.158180463556]
ϕ = electrostatic_potential(f, x, sim_box, rep_factors, sr_cutoff, k_rep_factors, α)
println("ϕ = ", ϕ * q_test)
println("should be -2676.8230141 K")
