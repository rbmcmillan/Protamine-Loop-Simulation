This is MATLAB software to accompany the article "DNA looping by protamine follows a nonuniform spatial distribution." 

sampleDist is a helper method that carries out the rejection sampling used in SimulateMultipleProtBinding_v6. This
file was authored by Dmitri Savransky. All attribution for this function should be given to him.

SimulateLoopFormation_v3gamma carries out the random looping model. See the comments at the top of the file for 
instructions on how to run the simulation.

SimulateLoopFormationWithElectroBias_v11 can simulate both the electrostatic binding and electrostatic multibinding
simulations. See comments at the top of the file for instructions on how to run the simulation. This function calls 
SimulateMultipleProtBinding_v6

SimulateMultipleProtBinding_v6 is a helper method that carries out one iteration of one or more protamines binding to 
a DNA molecule in an electrostatic potential.

DNAprop_withExtras_v3 is a MATLAB struct containing all experimental data. Loop start sites are saved in the field 
loop_start_site, while flower start sites are saved in the field flower_start_site.
