Base.@kwdef mutable struct Simulations_orig <: Simulations
    approach::String
    filename::String
    config_file::String
    SSPscen::Int64
    SSP::Int64
    nrun::Int64
    timehorizon::Int64
    intofuture::Int64
    totalimpact::Int64
    tmoney::Int64
    DR::Float64
    emissionyears::Array{Int64,1}
    FB_inv::Array{Float64,1}
    CO2_inv::Array{Float64,1}
    NOx_inv::Array{Float64,1}
    t0::Int64
    interp_years::Array{Int64}
    fueltype::Int64
    lc_analysis::Bool = false
    lc_CO2_inv::Array{Float64,1}
    lc_CH4_inv::Array{Float64,1}
    lc_N2O_inv::Array{Float64,1}
    lc_frac_CH4_ff::Array{Float64,1}
end

function getSimulationInfo_orig(input_dict::Dict)
    filename = input_dict["filename"];
    if !(typeof(filename) <: String)
        error("filename must be of type String")
    end

    config_file = input_dict["config_file"];
    if !(typeof(config_file) <: String)
        error("config_file must be of type String")
    end

    scenario = input_dict["scenario"];
    SSP_set = ["SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP4-3.4", "SSP4-6.0", "SSP5-3.4", "SSP5-8.5"]
    if !(scenario in SSP_set)
        error("scenario must be one of SSP scenarios in "*string(SSP_set))
    end
    SSPscen = findall(x->x == scenario, SSP_set)[1] #1-8, identifies which SSP-RCP scenario it is in full list of SSP_set
    SSP = SSPscenToSSP(SSPscen)                    #1-5, SSP number

    nrun = input_dict["nrun"];
    if !(typeof(nrun) <: Int)
        error("nrun must be an integer")
    end

    timehorizon = input_dict["timehorizon"];
    if !(typeof(timehorizon) <: Int)
        error("timehorizon must be an integer")
    end

    tmoney = input_dict["tmoney"];
    if !(typeof(tmoney) <: Int)
        error("tmoney must be an integer")
    elseif tmoney > 2021
        print("tmoney > 2021, changed to year 2021 USD")
        tmoney = 2021
    elseif tmoney < 1990
        print("tmoney < 1990, changed to year 1990 USD")
        tmoney = 1990
    end

    DR = input_dict["DR"];
    if DR < 0 || DR > 1
        error("DR must be in (0,1)")
    end

    emissionyears = input_dict["years"];
    if !(typeof(emissionyears) <: Array{Int,1})
        error("emissionyears must be an array of integers")
    end

    intofuture = length(emissionyears);
    totalimpact = timehorizon+intofuture;
    t0 = emissionyears[1];
    interp_years = interpYears(totalimpact, t0);

    fueltype = input_dict["fueltype"]
    if fueltype ∉ 1:16
        error("fueltype must be an integer between 1 and 16")
    end
    FB_inv = input_dict["kg_FB"]
    if !(typeof(FB_inv) <: Array{Float64,1}) || !(length(emissionyears) == length(FB_inv))
        error("kg_FB must be array of type Float with length = length(years)")
    end

    CO2_inv = input_dict["kg_CO2"]
    if !(typeof(CO2_inv) <: Array{Float64,1}) || !(length(emissionyears) == length(CO2_inv))
        error("kg_CO2 must be array of type Float with length = length(years)")
    end

    NOx_inv = input_dict["kg_NOx"]
    if !(typeof(NOx_inv) <: Array{Float64,1}) || !(length(emissionyears) == length(NOx_inv))
        error("kg_NOx must be array of type Float with length = length(years)")
    end

    lc_analysis = input_dict["lc_analysis"]
    if lc_analysis == 1
        lc_analysis = true
        lc_CO2_inv = input_dict["lc_kg_CO2"]
        lc_CH4_inv = input_dict["lc_kg_CH4"]
        lc_N2O_inv = input_dict["lc_kg_N2O"]
        lc_frac_CH4_ff = input_dict["lc_frac_CH4_ff"]
    elseif lc_analysis == 0
        lc_analysis = false
        lc_CO2_inv = zeros(length(emissionyears))
        lc_CH4_inv = zeros(length(emissionyears))
        lc_N2O_inv = zeros(length(emissionyears))
        lc_frac_CH4_ff = zeros(length(emissionyears))
    else
        error("lc_analysis must either be 0 (for false) or 1 (for true)")
    end

    simulation = Simulations_orig(
        approach="original",
        filename=filename,
        config_file=config_file,
        SSPscen=SSPscen,
        SSP=SSP,
        nrun=nrun,
        timehorizon=timehorizon,
        intofuture=intofuture,
        tmoney=tmoney,
        DR=DR,
        emissionyears=emissionyears,
        FB_inv=FB_inv,
        CO2_inv=CO2_inv,
        NOx_inv=NOx_inv,
        totalimpact=totalimpact,
        t0=t0,
        interp_years=interp_years,
        fueltype=fueltype,
        lc_analysis=lc_analysis,
        lc_CO2_inv=lc_CO2_inv,
        lc_N2O_inv=lc_N2O_inv,
        lc_CH4_inv=lc_CH4_inv,
        lc_frac_CH4_ff=lc_frac_CH4_ff,
        )
    return simulation
end