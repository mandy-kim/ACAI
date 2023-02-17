using Pkg
Pkg.activate(".")
Pkg.instantiate()

#Load ACAI
include("ACAI.jl")

#Specify input yaml file (make sure to check config_file and all other variables are correct)
input = "input/input_example_run_adv_approach.yml"

#Run ACAI
simulation, emission, climateoutput, aqoutput = runACAI(input);

#Outputs from run will be saved in output/$filename/, where the filename is specified in the input yml file
#If you rerun the above line (runACAI(input)), it will just read the outputs saved in this folder.
#If you want to perform a new run instead of reading previous outputs, change the (1) input or (2) filename in the input file, or just delete the output/$filename folder
