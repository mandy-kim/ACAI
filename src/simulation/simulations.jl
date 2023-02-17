abstract type Simulations end
include("simulations_adv.jl")
include("simulations_orig.jl")

function interpYears(totalimpact::Int64, t0::Int64)
    interp_years_array = collect(1:totalimpact) .+ (t0-1)
    return interp_years_array
end

function SSPscenToSSP(SSPscen::Int64)
    map = Dict{Int64,Int64}(
        1=>1,
        2=>1,
        3=>2,
        4=>3,
        5=>4,
        6=>4,
        7=>5,
        8=>5,
    )
    SSP = map[SSPscen]
    return SSP
end

"""
# getSimulationInfo - Simulation type constructor using user input 
### Input
- input_yaml_file: path to yaml file

### Outputs
- simulation
"""
function getSimulationInfo(input_yaml_file::String)
    input_dict = YAML.load_file(input_yaml_file)

    approach = input_dict["approach"];
    if approach ∉ ["advanced", "original"]
        error("approach must be specified as 'advanced' or 'original'")
    end

    if approach == "advanced"
        simulation = getSimulationInfo_adv(input_dict)
    else
        simulation = getSimulationInfo_orig(input_dict)
    end
    return simulation
end