"""
readACAIOutputs - read existing outputs for specific run/filename
"""
function readACAIOutputs(simulation::Simulations_adv)
    filename = simulation.filename
    f = jldopen("output/"*filename*"/outputs.jld2");
    emission = f["emission"]
    climateoutput = f["climateoutput"]
    aqoutput = f["aqoutput"]
    return emission, climateoutput, aqoutput
end

function readACAIOutputs(simulation::Simulations_orig)
    filename = simulation.filename
    f = jldopen("output/"*filename*"/outputs.jld2");
    emission = f["emission"]
    climateoutput = f["climateoutput"]
    return emission, climateoutput
end 

function runACAI_adv(simulation::Simulations_adv)
    println("running ACAI: advanced approach")
    #check if output folder already exists, just give outputs
    filename = simulation.filename;
    if isdir("output/"*"$filename/")
        println("reading existing output files")
        emission, climateoutput, aqoutput = readACAIOutputs(simulation)
        return emission, climateoutput, aqoutput
    end

    println("Initializing constants")
    #initialize variables uncertainty distributions
    climatedistribution, aqdistribution = MCDistributions(simulation);
    
    #read emissions and get annual sums in Emissions type (see emissionStructs.jl)
	emission = getYearlyEmis(simulation); #TOTAL annual emissions in kg and distance in km in 144 lon x 91 grid
    
    #initialize climate and health data constants
    climateconsts = getClimateConstants(simulation);
    healthconsts = getHealthConstants(simulation);
    
    #initialize struct to save outputs
    climateoutput =  ClimateOutputs(simulation);
    aqoutput = AirQualityOutputs(simulation, healthconsts);
    
    #read and save sensitivity data into variable (ONLY NOx FOR NOW, WILL NEED UPDATE ONCE OTHER SPECIES' DATA ARE READY)
    NOx_sens = getNOxSensitivities(simulation);
        #e.g., H2O_sens = ...
    
    #calculate NOx impacts (RF, AQ, Column O3) from emissions and sensitivities
    NOx_impacts, climateoutput, aqoutput = calculateNOxImpacts(NOx_sens, climateoutput, aqoutput, simulation, emission, healthconsts);
    
    #initialize contrail regions mask
    contrail_region_mask = load_region_map();

    println("Starting MC runs")
    for i in 1:simulation.nrun
        #Climate: RF -> ΔT
        #CO2 and Contrail RF
        climateoutput = update_output_RFs(i, climateoutput, simulation, emission, contrail_region_mask, climatedistribution, climateconsts);
        #CO2, Contrail, NOx, and background ΔT
        climateoutput = update_output_ΔTs(i, NOx_impacts, climateoutput, simulation, climatedistribution, climateconsts);
        
        #Air quality: ΔX -> ΔM -> $
        aqoutput = update_output_ΔMs(i, NOx_impacts, aqoutput, aqdistribution, simulation, healthconsts);
        aqoutput = update_output_ΔM_costs(i, aqoutput, aqdistribution, simulation, healthconsts);
    end
    println("MC runs done")
    climateoutput = update_output_ΔT_costs(climateoutput, simulation, climatedistribution, climateconsts);
    
    println("saving outputs in folder: "*"output/"*"$filename/")
    mkdir("output/"*"$filename/")
    jldsave("output/"*"$filename/outputs.jld2"; simulation, emission, climateoutput, aqoutput)
    
    return emission, climateoutput, aqoutput
end

function runACAI_orig(simulation::Simulations_orig)
    println("running ACAI: original approach")
    #check if output folder already exists, just give outputs
    filename = simulation.filename;
    if isdir("output/"*"$filename/")
        println("reading existing output files")
        emission, mainoutput = readACAIOutputs(simulation)
        return emission, mainoutput
    end

    println("Initializing constants")
    #initialize variables uncertainty distributions
    if simulation.lc_analysis 
        maindistribution, combustdistribution, lcdistribution = MCDistributions(simulation; lc_analysis=true);
    else
        maindistribution, combustdistribution = MCDistributions(simulation);
    end

    #get emission info for entire time horizon
    emission = getYearlyEmis(simulation);

    #initialize climate data constants
    climateconsts = getClimateConstants(simulation);
    
    #initialize struct to save outputs
	mainoutput = MainOutputs(simulation);
    if simulation.lc_analysis
        lcoutput = LCOutputs(simulation);
    end

    println("Starting MC runs")
    for i in 1:simulation.nrun
        #emissions -> RF -> ΔT
        mainoutput = update_output_RFs(i, mainoutput, simulation, emission, combustdistribution, maindistribution, climateconsts)
        mainoutput = update_output_ΔTs(i ,mainoutput, simulation, maindistribution, climateconsts)
        if simulation.lc_analysis
            lcoutput = update_output_RFs_lc(i, lcoutput, simulation, emission, lcdistribution, maindistribution, climateconsts)
            lcoutput = update_output_ΔTs_lc(i, lcoutput, simulation, maindistribution, climateconsts)
        end
    end
    println("MC runs done")
    #ΔT -> $
    mainoutput = update_output_ΔT_costs(mainoutput, simulation, maindistribution, climateconsts);
    if simulation.lc_analysis
        lcoutput = update_output_ΔT_costs_lc(lcoutput, mainoutput, simulation, maindistribution, climateconsts);
    end
    println("saving outputs in folder: "*"output/"*"$filename/")
    mkdir("output/"*"$filename/")
    if simulation.lc_analysis
        jldsave("output/"*"$filename/outputs.jld2"; simulation, emission, mainoutput, lcoutput)
        return emission, mainoutput, lcoutput
    else
        jldsave("output/"*"$filename/outputs.jld2"; simulation, emission, mainoutput)
        return emission, mainoutput
    end
end
