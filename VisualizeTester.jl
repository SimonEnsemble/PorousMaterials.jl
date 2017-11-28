using GLVisualize, GLWindow, PorousMaterials

filename = "ZASJAG_clean_min_charges.cif"

frame = read_crystal_structure_file(filename);

reps = replication_factors(frame.box, 14.0);

# Hard coding a He atom in
temp = zeros(3,1);
mol = Molecule(1,["He"], temp, [0.0]);

ljforcefield = read_forcefield_file("UFF.csv");

# Try to explore frame
#Energy,coord = exploreframe(frame, mol, ljforcefield, reps, mesh=41);
Energy = takesnapshot(frame, mol, ljforcefield, [10.,10.,10.] , reps, mesh=[21,21,21], startpoint=[15.,15.,15]);
#@printf("Minimum Energy = %f\nMaximum Energy = %f\n",minimum(Energy),maximum(Energy))
#for i=1:length(Energy)
#    if Energy[i] > 400
#        Energy[i] = 400
#    elseif Energy[i] < -400
#        Energy[i] = -400
#    end
#end

if !isdefined(:runtests)
    window = glscreen()
    timesignal = bounce(linspace(0f0,1f0,360))
end

description = """
Iso surface volume Plot.
"""


function volume_data(N)
    vol     = Energy+0
    max     = maximum(vol)
    min     = minimum(vol)
    vol     = (vol .- min) ./ (max .- min)
end

temp = Energy+0.0
vol = visualize(temp, :iso, isovalue=timesignal)
_view(vol, window)


if !isdefined(:runtests)
    renderloop(window)
end
