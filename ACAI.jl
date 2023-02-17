"""
ACAI: Aviation Climate and Air quality Impacts
"""
#Add packages required by ACAI
using Sobol
using Random
using YAML
using Distributions
using JLD2 
using NCDatasets
using Interpolations
using LinearAlgebra
using Statistics
using PyCall

#Load functions and modules to run ACAI
include("src/simulation/simulations.jl")
include("src/distribution/distributions.jl")
include("src/emission/emissions.jl")
include("src/constants/constants.jl")
include("src/output/outputs.jl")
include("src/sensitivities/sensitivities.jl")
include("src/climate/climate.jl")
include("src/airquality/airquality.jl")
include("src/main/ACAIfunctions.jl")

function runACAI(input::String)
    if !isdir("output/")
        mkdir("output/")
    end
	#read user inputs, save information in simulation variable/struct
	simulation = getSimulationInfo(input);
    #run appropriate model (advanced or original)
    if simulation.approach == "advanced" 
        emission, climateoutput, aqoutput = runACAI_adv(simulation)
        return simulation, emission, climateoutput, aqoutput
    elseif simulation.approach == "original"
        if simulation.lc_analysis
            emission, mainoutput, lcoutput = runACAI_orig(simulation)
            return simulation, emission, mainoutput, lcoutput
        else
            emission, mainoutput = runACAI_orig(simulation)
            return simulation, emission, mainoutput
        end
    end
end