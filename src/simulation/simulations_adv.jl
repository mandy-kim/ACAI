Base.@kwdef mutable struct Simulations_adv <: Simulations
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
    sensitivityfolders::Array{String}
    sensitivitycases::Int64
    sensitivityRFyears::Int64
    sensitivityAQyears::Int64
    sensitivityO3years::Int64
    emissionname::String
    emissionfolders::Array{String}
    emissionyears::Array{Int64,1}
    frequency::String
    weights::String
    layers::Int64
    t0::Int64
    interp_years::Array{Int64}
    VSLtype::String
    PM25_CRF::String
    O3_risk::String
    cessationlag::Bool
end

function getSimulationInfo_adv(input_dict::Dict)
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

    sensitivityfolders = input_dict["sensitivityfolders"];
    if !(typeof(sensitivityfolders) <: Array{String,1})
        error("sensitivityfolders must be of an array of strings")
    end

    sensitivitycases = input_dict["sensitivitycases"];
    if !(typeof(sensitivitycases) <: Int)
        error("sensitivitycases must be an integer")
    end

    sensitivityRFyears = input_dict["sensitivityRFyears"];
    if !(typeof(sensitivityRFyears) <: Int)
        error("sensitivityRFyears must be an integer")
    end

    sensitivityAQyears = input_dict["sensitivityAQyears"];
    if !(typeof(sensitivityAQyears) <: Int)
        error("sensitivityAQyears must be an integer")
    end

    sensitivityO3years = input_dict["sensitivityO3years"];
    if !(typeof(sensitivityO3years) <: Int)
        error("sensitivityO3years must be an integer")
    end

    emissionname = input_dict["emissionname"];
    if !(typeof(emissionname) <: String)
        error("emissionname must be a string")
    end

    emissionfolders = input_dict["emissionfolders"];
    if !(typeof(emissionfolders) <: Array{String,1})
        error("emissionfolders must be of an array of strings")
    end

    frequency = input_dict["frequency"];
    if frequency ∉ ["daily", "monthly", "yearly"]
        error("frequency must be 'daily', 'monthly', or 'yearly'")
    end

    layers = input_dict["layers"];
    if layers ∉ [47, 72]
        error("layers must be 47 or 72, corresponding to 47 or 72 vertical GEOS-chem grid of input emissions")
    end

    regrid_weights = input_dict["weights"];
    if ! isfile(regrid_weights)
        error("regrid_weights file does not exist")
    end

    VSLtype = input_dict["VSLtype"]
    if VSLtype ∉ ["country specific", "US", "global average"]
        error("VSLtype must be 'country specific', 'US', or 'global average'")
    end

    PM25_CRF = input_dict["PM25_CRF"]
    if PM25_CRF ∉ ["GEMM", "LL"]
        error("PM25_CRF must be either 'GEMM' or 'LL'")
    end

    O3_risk = input_dict["O3_risk"]
    if O3_risk ∉ ["respiratory", "all-cause"]
        error("O3_risk must be either 'respiratory' or 'all-cause'")
    end

    cessationlag = input_dict["cessationlag"]
    if cessationlag == 1
        cessationlag = true
    elseif cessationlag == 0
        cessationlag = false
    else
        error("cessationlag must either be 0 (for false) or 1 (for true)")
    end

    simulation = Simulations_adv(
        approach="advanced",
        filename=filename,
        config_file=config_file,
        SSPscen=SSPscen,
        SSP=SSP,
        nrun=nrun,
        timehorizon=timehorizon,
        intofuture=intofuture,
        tmoney=tmoney,
        DR=DR,
        sensitivityfolders=sensitivityfolders,
        sensitivitycases=sensitivitycases,
        sensitivityRFyears=sensitivityRFyears,
        sensitivityAQyears=sensitivityAQyears,
        sensitivityO3years=sensitivityO3years,
        emissionname=emissionname,
        emissionfolders=emissionfolders,
        emissionyears=emissionyears,
        frequency = frequency,
        weights=regrid_weights,
        layers=layers,
        totalimpact=totalimpact,
        t0=t0,
        interp_years=interp_years,
        VSLtype=VSLtype,
        PM25_CRF=PM25_CRF,
        O3_risk=O3_risk,
        cessationlag=cessationlag
        )
    return simulation
end